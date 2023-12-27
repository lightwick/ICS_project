%% Code Description
% (Modulation Diversity / Reduced Time) + (Spatial Modulation - Orthogonal Antenna Set)

% Real valued representation

%% Environment Variables
clear all
dbstop if error

% SM bits per iteration: bits sent per iteration
SM_bpi = 1;

TimeFrame = 2;

Nt = 8;
Nr = 8;
Np = 4;

EsN0_dB = 0:2:20;
EsN0 = db2pow(EsN0_dB);

eta = 2;

iteration = 1;

%% BEGIN
% eta/2 bits per dimension
% to represent k bits, we need 2^k levels
M = 2^(eta/2); % Modulation Order for one dimension

% EbN0 = db2pow(EbN0_dB);
% % EsN0 = EbN0 * bits per signal
% EsN0 = EbN0 * eta;
% EsN0_dB = pow2db(EsN0);

%     NormalizationFactor = sqrt(2/3*(M-1) * n);
NormalizationFactor = sqrt(2/3*(M^2-1)*Np);

if TimeFrame==2
    lambda = (1+sqrt(5))/2; % (1-sqrt(5))/2; both are plausible
    a = 1/sqrt(1+lambda^2);
    b = lambda*a;
    R = [a -b;
        b  a];
elseif TimeFrame==3
    phi = 2*cos(6/7*pi);
    alpha = (1+phi)/(1+phi+phi^2);
    beta = phi*alpha;
    gamma = -phi/(1+phi)*alpha;

    R = [alpha, beta, -gamma;
        beta, gamma, -alpha;
        gamma, alpha, -beta]
else
    error('TimeFrame has to be 2 or 3; others not implemented yet');
end

%% Candidate Generation
CandidateNum = M^(2*Np*TimeFrame); % Number of candidates
Candidates_MD = zeros(2*Np, TimeFrame, CandidateNum);
disp(sprintf("%d MD Candidates (before spatial modulation(SM) application)", size(Candidates_MD,3)));

for ii = 0:CandidateNum-1
    tmp = ii;
    count = 0;
    while tmp~=0
        Candidates_MD(floor(count/TimeFrame)+1, mod(count, TimeFrame)+1, ii+1) = mod(tmp, M);
        tmp = floor(tmp/M);
        count = count + 1;
    end
end

CandidateSymbol = pammod(Candidates_MD, M) / NormalizationFactor;

CandidateSymbol = pagetranspose(pagemtimes(R,'none', CandidateSymbol, 'transpose'));

%% Candidate Gen, restart
for page=1:CandidateNum
    for ii=2:TimeFrame
        CandidateSymbol(:, ii, page) = circshift(CandidateSymbol(:, ii, page), ii-1);
    end
end
%     CandidateSymbol = pagetranspose(CandidateSymbol);
% Resize Page
% CandidateSymbol = reshape(CandidateSymbol, [], 1, CandidateNum);

BEC = zeros(length(EsN0), 1);
EuclideanDistance = zeros(CandidateNum, 1);

%% Validating Normalization
TotalSum = zeros(size(CandidateSymbol,1), size(CandidateSymbol, 2));

for page=1:CandidateNum
    TotalSum(:,:) = TotalSum(:,:) + (CandidateSymbol(:,:,page).^2);
end


%% Spatial Modulation (Orthogonal Antenna Set)
codebook = getCodebook();
% no rotation for now

assert(Nt==2*Np, "The current implementation requires Nt=2*Np");
wordperbook = 2; % Because Nt/Np=2 for in this implementation

% HARD CODED: there should be a better way to implement this.
if SM_bpi==1
    SM_Candidates = [1 2; 2 1];
elseif SM_bpi==2
    SM_Candidates = [1 2; 2 1; 3 4; 4 3];
else
    error(['Only SM_bpcu of 1 or 2 implemented.\n' ...
        'If need more, change codebook']);
end

% c is the number of antenna combinations when applying SM
c = length(SM_Candidates);
Candidates_SM = zeros(Nt*2, size(Candidates_MD, 2), length(Candidates_MD)*c);

TransmitCandidate = zeros(length(Candidates_SM), TimeFrame);

count = 1;
for i1 = 1 : length(SM_Candidates)
    [BookIndex, WordIndex] = transmit2index(SM_Candidates(i1, :), wordperbook);
    for i3=1:length(CandidateSymbol)
        for TimeSlot = 1:TimeFrame
            ApplyingAntenna = codebook{BookIndex(TimeSlot)}{WordIndex(TimeSlot)};
            for ChoiceIndex=1:Np
                Candidates_SM([ApplyingAntenna(ChoiceIndex), ApplyingAntenna(ChoiceIndex)+Nt], TimeSlot, count) = CandidateSymbol([ChoiceIndex, ChoiceIndex+Np], TimeSlot, i3);
            end
        end
        TransmitCandidate(count, :) = SM_Candidates(i1, :) - 1;
        count = count + 1; % just simply lazy coding
    end
end

%% Main simulation process
% Timer
FivePercent = ceil(iteration/20);

for iTotal = 1:iteration
    if mod(iTotal-100, FivePercent)==0
        tic
    end
    %% Bit Generation
    % Transmitter Defined Bits
    TransmitterDecimal = randi([0 c-1], 1, 1); % Two Time Frame
    TransmitterBinary = de2bi(TransmitterDecimal, log2(c), 'left-msb');

    AntChoice = SM_Candidates(TransmitterDecimal+1, :);

    [CodebookIndex, CodewordIndex] = transmit2index(AntChoice, wordperbook);
    % CodebookIndex = floor(TransmitterDecimal/a)+1;
    % CodewordIndex = mod(TransmitterDecimal, a)+1;

    SignalSequence = randi([0 M-1], 2*Np, TimeFrame);
    % SignalBinary = de2bi(SignalSequence, log2(M), 'left-msb');
    % SymbolSequence = qammod(SignalSequence, M) / NormalizationFactor;
    SymbolSequence = pammod(SignalSequence, M) / NormalizationFactor;

    RotatedSymbol = (R*SymbolSequence')';
    for ii=2:TimeFrame
        RotatedSymbol(:, ii) = circshift(RotatedSymbol(:,ii), ii-1);
    end
    % x2_r = RotatedSymbol(:);

    RotatedSymbol_SM = zeros(Nt*2, TimeFrame);

    %% Applying Spatial Modulation
    for TimeSlot=1:TimeFrame
        AntennaChoice = codebook{CodebookIndex(TimeSlot)}{CodewordIndex(TimeSlot)};
        for ChoiceIndex=1:Np
            ApplyingAntenna = AntennaChoice(ChoiceIndex);
            RotatedSymbol_SM([ApplyingAntenna, ApplyingAntenna+Nt], TimeSlot) = RotatedSymbol([ChoiceIndex, ChoiceIndex+Np], TimeSlot);
        end
    end

    % x2_r = RotatedSymbol_SM(:); % Real valued representation
    x2_r = RotatedSymbol_SM;
    
    %% Coding
    H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) / sqrt(2); % Receiver x Transmitter
    H_r = [real(H), -imag(H);
                  imag(H), real(H)];

    Noise = randn(size(H_r, 1), 2) / sqrt(2);

    %% DEBUG
    % Noise = zeros(size(Noise));

    Candidate_y = pagemtimes(H_r, Candidates_SM);

    for SNR_idx = 1 : length(EsN0)
        ReceivedSignal = H_r * x2_r + Noise / sqrt(EsN0(SNR_idx));
        
        % Decoding
        EuclideanDistance = pagenorm(Candidate_y - ReceivedSignal, 'fro').^2;
        [~, idx] = min(EuclideanDistance);
        
        DetectedTransmitter = TransmitCandidate(idx, 1);
        DetectedSignal = Candidates_MD(:, :, mod(idx-1, length(CandidateSymbol))+1);
        DetectedTransmitterBinary = de2bi(DetectedTransmitter, log2(c), 'left-msb');

         % NOTE: THIS ONLY WORKS BECAUSE THE MODULATION ORDER IN PAMMOD IS 2; MEANING ONLY 0 AND 1 IS INSIDE THE 'DetectedTransmitter' and 'DetectedSignal' variable.
        TransmitError = sum((TransmitterBinary~=DetectedTransmitterBinary), 'all');
        SignalError = sum((de2bi(SignalSequence, log2(M))~=de2bi(DetectedSignal, log2(M))), 'all');
        ErrorCount = TransmitError + SignalError;
        
        BEC(SNR_idx) = BEC(SNR_idx) + ErrorCount;
        % TBEC(1, SNR_idx) = TBEC(1, SNR_idx) + TransmitError;
        % SBEC(1, SNR_idx) = SBEC(1, SNR_idx) + SignalError;
        
        % if ErrorCount~=0
        %     disp(sum((SignalBinary~=DetectedBinary), 'all'))
        %     disp(sum((TransmitterBinary~=DetectedTransmitterBinary), 'all'))
        % end
    end

    if mod(iTotal-100, FivePercent)==0
        ElapsedTime = toc;
        EstimatedTime = (iteration-iTotal)*ElapsedTime;
        disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/iteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
    end
end

%% Need to fix 분모 when calculating BER
BitsPerIteration = TimeFrame * (log2(M)*Np*2) + SM_bpi;
TotalTransmitBits = BitsPerIteration * iteration;

BER = (BEC/TotalTransmitBits)';
% TBER(1,:) = TBEC(1,:)/(log2(c)*iteration);
% SBER(1,:) = SBEC(1,:)/(2*log2(M)*iteration);

%% Plotting
% BER_Title = sprintf("BER for %d-QAM", M);
BER_Title = "Orthogonal MD-SM: N_T=8, N_p=4"
x_axis = "E_s/N_0 (dB)";

legend_order = [""];
myplot(EsN0_dB, BER, BER_Title, x_axis, "BER", legend_order);
ylim([10^(-6) 1])


%% Helper functions
% refer to codebook_gen.m file
function codebook = getCodebook()
    codebook = {
        {[1,2,3,4], [5,6,7,8]};
        {[1,2,7,8], [3,4,5,6]}
        };
end

function [BookIndex, WordIndex] = transmit2index(index, wordperbook)
    BookIndex = ceil(index/wordperbook);
    WordIndex = mod(index-1, wordperbook)+1;
end