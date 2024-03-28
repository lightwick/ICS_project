%close all
clear all
clc

%% Code Description
% Regular Modulation Diversity + Direct Extension SM
% Addition to MD_SM_DirectExtension to work on R3

%% TODO
% check if input configuration is viable with the code
% condition 1: in M=16, n_t<=8

% 24/03/28 Remove 'left-msb' in de2bi functions; really redundant

%% Current STATUS
% When Nt=2, no BER==0; meaning there is a problem when SM is added

%% DEBUGGING FEATURES
DEBUG = false;

%% User Defined Configuration
addpath('../../tools');
M = 2;
EsN0_dB =  0:2:20;
EsN0 = db2pow(EsN0_dB);

Nt = 4;
Nr = Nt;
TimeFrame = 3;
Np = TimeFrame; % Np being the number of antennas that transmit signals; another thought

iteration = 10^4;

assert(Nt==4,'Only implemented Nt=4');

%% Basic Configuration

num = 1;
BEC = zeros(num, length(EsN0_dB));
TBEC = zeros(num, length(EsN0_dB));
SBEC = zeros(num, length(EsN0_dB));

NormalizationFactor = sqrt(2/3*(M^2-1) * Np); % Normally, this would be 'NormalizationFactor = sqrt(2/3*(M^2-1) * Nt);', but Spatial Modulation is added here


%% Modulation Diversity; Setup
if TimeFrame==2
    lambda = (1+sqrt(5))/2;
    % lambda = (1-sqrt(5))/2;
    a = 1/sqrt(1+lambda^2);
    b = lambda*a;
    R = [a -b;
           b  a];
elseif TimeFrame==3
    phi = 2*cos(6*pi/7);
    alpha = (1+phi)/(1+phi+phi^2);
    beta = phi*alpha;
    gamma = -phi/(1+phi)*alpha;
    R = [alpha beta -gamma; beta gamma -alpha; gamma alpha -beta];
else
    error('TimeFrame must be 2 or 3');
end

% Candidate Generation
CandidateNum = M^(2*Np^2) % * nchoosek(Nt, 2); % Number of candidates
BeforeModulation = zeros(2*Np, Np, CandidateNum);
disp(sprintf("%d MD Candidates (before spatial modulation(SM) application)", size(BeforeModulation,3)));

for ii = 0:CandidateNum-1
    TransmitterNum = ii;
    count = 0;
    while TransmitterNum~=0
        BeforeModulation(floor(count/TimeFrame)+1, mod(count, TimeFrame)+1, ii+1) = mod(TransmitterNum, M);
        TransmitterNum = floor(TransmitterNum/M);
        count = count + 1;
    end
end

CandidateSymbol = pammod(BeforeModulation, M) / NormalizationFactor;
CandidateSymbol = pagetranspose(pagemtimes(R,'none', CandidateSymbol, 'transpose'));

%% Candidate Gen, restart
for page=1:CandidateNum
    for ii=2:TimeFrame
        CandidateSymbol(:, ii, page) = circshift(CandidateSymbol(:,ii, page), ii-1);
    end
end

%% Spatial Modulation; Setup
% map = codebook_gen();
% codebook=map(Nt);

codebook = nchoosek([1:Nt], Np);
c = length(codebook);

if Nt<=4
    switch M
        case 2
            theta(2) = 1.57
        case 4
            theta(2) = 0.61
        case 16
            theta(2) = 9.05
        case 64
            theta(2) = 8.23
        otherwise
            error('M value needs to be among the following: 2, 4, 16, 64')
    end
else
    for k=2:n
        switch M
            case 2
                theta(k)=(k-1)*pi/(n);
            case {4, 16}
                theta(k)=(k-1)*pi/(2*n);
            otherwise
                error('M value needs to be 2,4, or 16 in n_t>4 case')
        end
    end
end

%% Creating SM Candidates
SMAppliedCandidate = zeros([size(CandidateSymbol) length(codebook)]);

%% ????? 이거 되는 게 맞는지 확인 해야함;;;;
for i1=1:length(codebook)
    SMAppliedCandidate([codebook(i1, :) codebook(i1, :)+Nt], :, :, i1) = CandidateSymbol(:, :, :); %% 이게 돼????
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

    % AntChoice = SM_Candidates(TransmitterDecimal+1, :);
    AntChoice = codebook(TransmitterDecimal+1, :);

    % [CodebookIndex, CodewordIndex] = transmit2index(AntChoice, a);

    SignalSequence = randi([0 M-1], 2*Np, Np);

    SymbolSequence = pammod(SignalSequence, M) / NormalizationFactor;

    RotatedSymbol = (R*SymbolSequence')';
    for ii=2:Np
        RotatedSymbol(:, ii) = circshift(RotatedSymbol(:,ii), ii-1);
    end
    % x2_r = RotatedSymbol(:);

    RotatedSymbol_SM = zeros(Nt*2, Np);
    RotatedSymbol_SM([AntChoice AntChoice+Nt], :) = RotatedSymbol; %% 이게 돼????

    x2_r = RotatedSymbol_SM;
    
    %% Coding
    H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) / sqrt(2); % Receiver x Transmitter
    H_r = [real(H), -imag(H);
                  imag(H), real(H)];

    Noise_r = randn(2*Nr, TimeFrame) / sqrt(2);

    % DEBUG
    Noise_r = zeros(size(Noise_r));

    % Candidate_y = pagemtimes(H_r, Candidates_SM);
    Candidate_y = pagemtimes(H_r, SMAppliedCandidate);

    for SNR_idx = 1 : length(EsN0)
        ReceivedSignal = H_r * x2_r + Noise_r / sqrt(EsN0(SNR_idx));
        
        % Decoding
        EuclideanDistance = pagenorm(Candidate_y - ReceivedSignal, 'fro').^2;
        [val, idx] = min(EuclideanDistance);
        [~, TransmitterNum] = min(val);

        DetectedSM = TransmitterNum-1;
        DetectedSignal = BeforeModulation(:,:,idx(TransmitterNum));
        DetectedTransmitterBinary = de2bi(DetectedSM, log2(c), 'left-msb');

         % NOTE: THIS ONLY WORKS BECAUSE THE MODULATION ORDER IN PAMMOD IS 2; MEANING ONLY 0 AND 1 IS INSIDE THE 'DetectedTransmitter' and 'DetectedSignal' variable.
        TransmitError = sum((TransmitterBinary~=DetectedTransmitterBinary), 'all');
        SignalError = sum((de2bi(SignalSequence, log2(M))~=de2bi(DetectedSignal,log2(M))), 'all');
        ErrorCount = TransmitError + SignalError;
        

        BEC(1, SNR_idx) = BEC(1, SNR_idx) + ErrorCount; %bitcount(bitxor([detected_x1; detected_x2], SignalSequence)));
        TBEC(1, SNR_idx) = TBEC(1, SNR_idx) + TransmitError; %bitcount(bitxor([detected_x1; detected_x2], SignalSequence)));
        SBEC(1, SNR_idx) = SBEC(1, SNR_idx) + SignalError; %bitcount(bitxor([detected_x1; detected_x2], SignalSequence)));
        
        % if ErrorCount~=0
        %     disp(sum((SignalBinary~=DetectedBinary), 'all'))
        %     disp(sum((TransmitterBinary~=DetectedTransmitterBinary), 'all'))
        % end
    end

    if mod(iTotal-100, FivePercent)==0
        ElapsedTime = toc;
        EstimatedTime = (iteration-iTotal)*ElapsedTime;
        hours = floor(EstimatedTime / 3600);
        minutes = floor(mod(EstimatedTime, 3600) / 60);
        seconds = floor(mod(EstimatedTime, 60));
        disp(sprintf("%d%%, estimated wait time %d hours %d minutes %d seconds", round(iTotal/iteration*100), hours, minutes, seconds))
    end
end

%% Need to fix 분모 when calculating BER
%% 돌아와서 다시 확인
BitsPerIteration = TimeFrame * (log2(M)*Np*2) + log2(c);
TotalTransmitBits = BitsPerIteration * iteration;

BER(1,:) = BEC(1,:)/TotalTransmitBits;
% TBER(1,:) = TBEC(1,:)/(log2(c)*iteration);
% SBER(1,:) = SBEC(1,:)/(2*log2(M)*iteration);

%% Plotting
% BER_Title = sprintf("BER for %d-QAM", M);
BER_Title = "Orthogonal MD-SM: N_T=4, N_p=2"
x_axis = "E_s/N_0 (dB)";

legend_order = ["4-QAM Equivalent Spectral Efficiency; 5 bpcu", "16-QAM Equivalent Spectral Efficiency; 9 bpcu", "MD-SM Reduced Time Slot"];
myplot(EsN0_dB, BER, BER_Title, x_axis, "BER", legend_order);
ylim([10^(-6) 1])


%% Tools
function Candidates = get_candidates(M, Nt)
    AllNumbers = de2bi([0:M^Nt-1], Nt*log2(M), 'left-msb');
    Candidates = zeros(M^Nt, Nt);
    for ii = 1 : M^Nt
        for jj = 1 : Nt
            Candidates(ii,jj) = bi2de(AllNumbers(ii,log2(M)*(jj-1)+1:log2(M)*jj), 'left-msb');
        end
    end
end

function [BookIndex, WordIndex] = transmit2index(index, wordperbook)
    BookIndex = ceil(index/wordperbook);
    WordIndex = mod(index-1, wordperbook)+1;
end