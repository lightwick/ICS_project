close all
clear all
clc

%% TODO

%% Current STATUS
% When Nt=2, no BER==0; meaning there is a problem when SM is added

%% User Defined Configuration
addpath('../../tools');

% NONOISE = input("No noise?");

EsN0_dB =  0:2:20;
% EsN0_dB = 0;
EsN0 = db2pow(EsN0_dB);

Nt = 4;
Np = 2;
Nr = 4;
M = 4;

iteration = 10^4;

num = 1;
BEC = zeros(num, length(EsN0_dB));
TBEC = zeros(num, length(EsN0_dB));
SBEC = zeros(num, length(EsN0_dB));

NormalizationFactor = sqrt(2/3*(M-1) * Np); % Normally, this would be 'NormalizationFactor = sqrt(2/3*(M-1) * Nt);', but Spatial Modulation is added here

% map = codebook_gen();
% codebook=map(Nt);

c = Nt*(Nt-1)/2;
while bitand(c, c-1) % a way to check if a number is the power of 2
    c=c-1;
end

assert(Nt==4 && Np==2, "Nt=4, Np=2 only");
codebook = [[1,2]; [3,4]; [1,3]; [2,4]];

%% Modulated Candidate Generation
% Candidates = qammod(0:M-1, M) / NormalizationFactor;
Candidates = get_candidates(M, Np);
ModulatedCandidates = qammod(Candidates, M) / NormalizationFactor;
assert(M==4 || M==16, 'Constraint on using current normalization factor for candidate generation');

%% Main simulation process
% Timer
FivePercent = ceil(iteration/20);

for iTotal = 1:iteration
    if mod(iTotal-100, FivePercent)==0
        tic
    end
    %% Bit Generation
    % Transmitter Defined Bits
    TransmittedCodeword_idx = randi([1 c]); % 1-based index
    AntSelection = codebook(TransmittedCodeword_idx, :);

    H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) / sqrt(2);
    n = (randn(Nr, 1) + 1j * randn(Nr, 1)) / sqrt(2);
    SignalDecimal = randi([0 M-1], Np, 1);

    theta = exp(1j*[0, pi/8, pi/4, 3*pi/8]);
    x = qammod(SignalDecimal, M, 'UnitAveragePower', true) / sqrt(Np) * theta(TransmittedCodeword_idx);
    
    x_eff = zeros(Nt, 1);
    x_eff(AntSelection) = x;
    %% SM Application

    % n = zeros(size(n));

    for SNR_idx = 1 : length(EsN0)
        y = H * x_eff + n / sqrt(EsN0(SNR_idx));
        foo = inv(H' * H) * H';
        [val, idx] = sort(abs(foo*y), 'descend');

        % for i1=1:length(codebook)
        %     T_k(:,:,i1) = H(:, codebook(i1, :));
        % end
        % T_k_hermitian = pagectranspose(T_k);
        % foo = pageinv(pagemtimes(T_k_hermitian, T_k));
        % foo = pagemtimes(foo, T_k_hermitian);        
        % foo = pagemtimes(foo, y);
        
        % [val, idx] = sort(abs(y-foo), 'descend');
        idx = sort(idx(1:Np));

        for i1 = 1:c
            if codebook(i1,:)==idx.'
                DetectedWord_idx = i1;
                break
            end
        end

        Candidate_y = zeros(Nt, length(ModulatedCandidates));
        for i1=1:length(ModulatedCandidates)
            Candidate_y(idx, i1) = ModulatedCandidates(i1,:);
        end
        Candidate_y = H * Candidate_y * theta(DetectedWord_idx);
        diff = sum(abs(y-Candidate_y).^2, 1);
        [val, idx] = min(diff);

        DetectedBinary = de2bi(Candidates(idx, :), log2(M));
        SignalBinary = de2bi(SignalDecimal, log2(M));
        SignalError = sum((SignalBinary~=DetectedBinary), 'all');
        
        TransmitError = sum((de2bi(TransmittedCodeword_idx-1,log2(c))~=de2bi(DetectedWord_idx-1, log2(c))), 'all');
        ErrorCount = TransmitError + SignalError;
        
        BEC(1, SNR_idx) = BEC(1, SNR_idx) + ErrorCount; %bitcount(bitxor([detected_x1; detected_x2], SignalSequence)));
        TBEC(1, SNR_idx) = TBEC(1, SNR_idx) + TransmitError; %bitcount(bitxor([detected_x1; detected_x2], SignalSequence)));
        SBEC(1, SNR_idx) = SBEC(1, SNR_idx) + SignalError; %bitcount(bitxor([detected_x1; detected_x2], SignalSequence)));
    end

    if mod(iTotal-100, FivePercent)==0
        ElapsedTime = toc;
        EstimatedTime = (iteration-iTotal)*ElapsedTime;
        disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/iteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
    end
end

TotalTransmitBits = (Np * log2(M) + log2(c))* iteration;

BER(1,:) = BEC(1,:)/TotalTransmitBits;
TBER(1,:) = TBEC(1,:)/(log2(c)*iteration);
SBER(1,:) = SBEC(1,:)/(2*log2(M)*iteration);

%% Plotting
BER_Title = "STBC-SM";
x_axis = "SNR (dB)";

proposed_legend = sprintf("STBC-SM, n_T=%d, %d-QAM", Nt, M);
legend_order = ["STBC-SM 64-QAM(7 bpcu)", "MD-SM Orthogonal(9 bpcu)", "STBC-SM 16-QAM(5 bpcu)"];
myplot(EsN0_dB, BER, BER_Title, x_axis, "BER", legend_order);
ylim([10^(-6) 1])