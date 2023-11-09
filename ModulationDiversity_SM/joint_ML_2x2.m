%close all
clear all
clc

%% Verdict
% Joint search over 2 time slots is meaningless.
% There is no correlation

%% TODO
% check if input configuration is viable with the code
% condition 1: in M=16, n_t<=8

%% Current STATUS
% When Nt=2, no BER==0; meaning there is a problem when SM is added

%% DEBUGGING FEATURES
DEBUG = false;

%% User Defined Configuration
addpath('../tools');
M = 4;
EsN0_dB =  0:2:24;
EsN0 = db2pow(EsN0_dB);

Nt = 2;
Nr = 2;

iteration = 10^6;

%% Basic Configuration
num = 1;
BEC = zeros(num, length(EsN0_dB));

TimeFrame = 2;

% Candidate Generation
CandidateNum = M^(Nt*2);
Candidates = zeros(Nt, TimeFrame, CandidateNum);
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

CandidateSymbol = qammod(Candidates, M, 'UnitAveragePower', true) / sqrt(Nt);

%% Main simulation process
% Timer
FivePercent = ceil(iteration/20);

for iTotal = 1:iteration
    if mod(iTotal-100, FivePercent)==0
        tic
    end
    %% Bit Generation
    SignalSequence = randi([0 M-1], Nt, TimeFrame);
    SignalBinary = de2bi(SignalSequence(:), log2(M), 'left-msb');
    SymbolSequence = qammod(SignalSequence, M, 'UnitAveragePower', true) / sqrt(Nt);

    %% Coding
    H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) / sqrt(2); % Receiver x Transmitter
    Noise = (randn(Nr, TimeFrame) + 1j*randn(Nr, TimeFrame)) / sqrt(2);

    Candidate_y = pagemtimes(H, CandidateSymbol);

    for SNR_idx = 1 : length(EsN0)
        ReceivedSignal = H * SymbolSequence + Noise / sqrt(EsN0(SNR_idx));
        
        % Decoding
        EuclideanDistance = pagenorm(Candidate_y - ReceivedSignal, 'fro').^2;
        [~, idx] = min(EuclideanDistance);
        
        DetectedSignal = Candidates(:, :, idx);
        DetectedBinary = de2bi(DetectedSignal(:), log2(M), 'left-msb');

        ErrorCount = sum(DetectedBinary~=SignalBinary, 'all');

        BEC(1, SNR_idx) = BEC(1, SNR_idx) + ErrorCount; %bitcount(bitxor([detected_x1; detected_x2], SignalSequence)));
    end

    if mod(iTotal-100, FivePercent)==0
        ElapsedTime = toc;
        EstimatedTime = (iteration-iTotal)*ElapsedTime;
        disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/iteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
    end
end

%% Need to fix 분모 when calculating BER
%% 돌아와서 다시 확인
BitsPerIteration = Nt*log2(M);

BER = BEC / (BitsPerIteration*iteration);

%% Plotting
BER_Title = sprintf("BER for %d-QAM", M);
x_axis = "SNR (dB)";

legend_order = ["2x2 joint decoding"];
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