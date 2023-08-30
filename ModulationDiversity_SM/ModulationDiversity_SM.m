%close all
clear all
clc

%% TODO
% check if input configuration is viable with the code
% condition 1: in M=16, n_t<=8

%% Current STATUS
% When Nt=2, no BER==0; meaning there is a problem when SM is added

%% DEBUGGING FEATURES
DEBUG = false;
NOISE = true;

%% User Defined Configuration
addpath('../tools');
M = 2;
EsN0_dB =  0:2:22;
if DEBUG==true
    EsN0_dB=0;
end

Nt = 4;
Nr = Nt;
Np = 2; % Np being the number of antennas that transmit signals; another thought;

iteration = 10^5;

assert(Nt==4,'Only implemented Nt=4');

%% Basic Configuration
EsN0_config % configuration script in tools

num = 1;
BEC = zeros(num, length(EsN0_dB));
TBEC = zeros(num, length(EsN0_dB));
SBEC = zeros(num, length(EsN0_dB));

NormalizationFactor = sqrt(2/3*(M-1) * 2); % Normally, this would be 'NormalizationFactor = sqrt(2/3*(M-1) * Nt);', but Spatial Modulation is added here


%% Modulation Diversity; Setup
lambda = (1+sqrt(5))/2;
% lambda = (1-sqrt(5))/2;
a = 1/sqrt(1+lambda^2);
b = lambda*a;
R = [a -b;
       b  a];

TimeFrame = 2;

% Candidate Generation
CandidateNum = M^(2*Np^2) % * nchoosek(Nt, 2); % Number of candidates
Candidates = zeros(2*Np, Np, CandidateNum);
disp(sprintf("%d Candidates", size(Candidates,3)));

for ii = 0:CandidateNum-1
    tmp = ii;
    count = 0;
    while tmp~=0
        Candidates(floor(count/TimeFrame)+1, mod(count, TimeFrame)+1, ii+1) = mod(tmp, M);
        tmp = floor(tmp/M);
        count = count + 1;
    end
end

CandidateSymbol = pammod(Candidates, M) / NormalizationFactor;
CandidateSymbol = pagetranspose(pagemtimes(R,'none', CandidateSymbol, 'transpose'));

%% Candidate Gen, restart
for page=1:CandidateNum
    for ii=2:TimeFrame
        CandidateSymbol(:, ii, page) = circshift(CandidateSymbol(:,ii, page), ii-1);
    end
end
%     CandidateSymbol = pagetranspose(CandidateSymbol);
% Resize Page
% CandidateSymbol = reshape(CandidateSymbol, [], 1, CandidateNum);

%% Spatial Modulation; Setup
map = codebook_gen();
codebook=map(Nt);

c = Nt*(Nt-1)/2;
while bitand(c, c-1) % a way to check if a number is the power of 2
    c=c-1;
end 
% number of total codewords
n = size(codebook,1); % number of codebooks
a = size(codebook{1},2); % number of codewords per codebook

last_a = c-a*(n-1); % Has no purpose in this version of implementation(n_t=3,4,8)

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

%% Modulated Candidate Generation
% Candidates = zeros(M^2, 2);
% 
%  for ii = 0:M^2-1
%      for jj = 1:2
%          Candidates(ii+1,jj) = mod(floor(ii/M^(2-jj)),M);
%      end
%  end
%  Candidates = qammod(Candidates', M) / NormalizationFactor; % Each column is a Candidate

%% Creating SM Candidates
Candidates_SM = zeros(size(Candidates, 1), size(Candidates, 2), length(Candidates)*c^2);

count = 1;
for i1=1:n  % n is codebook
    for i2=1:a % a is word per book
        % I want to cry.....
        % My life has no meaning....
        % Fuck Everything
        % I want to die
        % 그래 슈밤바 그냥 이번 주말에 뭘 할지나 고민하자..
        % 슈밤 뭐하지.........

        %% NEED TO DO FROM HERE
        for TimeSlot=1:Np
            AntennaChoice = codebook{CodebookIndex(TimeSlot)}{CodewordIndex(TimeSlot)};
            for ChoiceIndex=1:Np
                ApplyingAntenna = AntennaChoice(ChoiceIndex);
                RotatedSymbol_SM([ApplyingAntenna, ApplyingAntenna+Nt], TimeSlot) = RotatedSymbol([ChoiceIndex, ChoiceIndex+Np], TimeSlot);
            end
        end
        count = count+1; % just fucking lazy implementation
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
    TransmitterDecimal = randi([0 c-1], 1, 2); % Two Time Frame
    TransmitterBinary = de2bi(TransmitterDecimal, log2(c), 'left-msb');

    CodebookIndex = floor(TransmitterDecimal/a)+1;
    CodewordIndex = mod(TransmitterDecimal, a)+1;

    SignalSequence = randi([0 M-1], 2*Np, Np);
    % SignalBinary = de2bi(SignalSequence, log2(M), 'left-msb');
    % SymbolSequence = qammod(SignalSequence, M) / NormalizationFactor;
    SymbolSequence = pammod(SignalSequence, M) / NormalizationFactor;

    RotatedSymbol = (R*SymbolSequence')';
    for ii=2:Np
        RotatedSymbol(:, ii) = circshift(RotatedSymbol(:,ii), ii-1);
    end
    % x2_r = RotatedSymbol(:);

    RotatedSymbol_SM = zeros(Nt*2, Np);

    %% Applying Spatial Modulation
    for TimeSlot=1:Np
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

    for SNR_idx = 1 : length(EsN0)
        ReceivedSignal = H_r * x2_r + Noise / sqrt(EsN0(SNR_idx));
        
        % Decoding


        TransmitError = sum((TransmitterBinary~=DetectedTransmitterBinary), 'all');
        SignalError = sum((SignalBinary~=DetectedBinary), 'all');
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
        disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/iteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
    end
end
TotalTransmitBits = (2 * log2(M) + log2(c))* iteration;

BER(1,:) = BEC(1,:)/TotalTransmitBits;
TBER(1,:) = TBEC(1,:)/(log2(c)*iteration);
SBER(1,:) = SBEC(1,:)/(2*log2(M)*iteration);

%% Plotting
BER_Title = sprintf("BER for %d-QAM", M);
x_axis = "SNR (dB)";

legend_order = ["STBC-SM, n_T=8, 16-QAM, 6 bits/s/Hz"];
myplot(EsN0_dB, BER, BER_Title, x_axis, "BER", legend_order);
ylim([10^(-6) 1])