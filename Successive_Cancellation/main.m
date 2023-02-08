close all
clear all
clc
dbstop if error
dbstop if warning

addpath('../tools/')

% Environment Varible
M = 16
Nt = 2
Nr = 4
NumberIteration = 10^4;

% Simulation
LengthBitSequence = Nt * log2(M); % log2(M) bits per signal
LengthSignalSequence = Nt;

% EbN0_dB = 0:5:25;
% EbN0 = db2pow(EbN0_dB);
% 
% EsN0 = EbN0 * log2(M);
% EsN0_db = pow2db(EsN0);

EsN0_dB = 0:5:25;
EsN0 = db2pow(EsN0_dB);

EbN0 = EsN0 / log2(M);
EbN0_dB = pow2db(EbN0);

% BitErrorCount_ZF = zeros(1, length(EsN0_dB));
% SignalErrorCount_ZF = zeros(1, length(EsN0_dB));
% BitErrorCount_MLD = zeros(1, length(EsN0_dB));
% SignalErrorCount_MLD = zeros(1, length(EsN0_dB));
% 
% BitErrorCount_MMSE = zeros(1, length(EsN0_dB));
% SignalErrorCount_MMSE = zeros(1, length(EsN0_dB));

BitErrorCount = zeros(7, length(EsN0_dB));
SignalErrorCount = zeros(7, length(EsN0_dB));

NormalizationFactor = sqrt(2/3*(M-1)*Nt);

FivePercent = ceil(NumberIteration/20);
for iTotal = 1 : NumberIteration
    if mod(iTotal-100, FivePercent)==0
        tic
    end
    % Bit Generation
    SignalSequence = randi([0 M-1], Nt, 1);
    SignalBinary = de2bi(SignalSequence, log2(M), 'left-msb');
    SymbolSequence = qammod(SignalSequence, M) / NormalizationFactor;
    
    NoiseSequence = (randn(Nr, 1) + 1j * randn(Nr, 1)) / sqrt(2); % Noise (n) Generation
    H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) ./ sqrt(2); % Receiver x Transmitter
    for indx_EbN0 = 1 : length(EsN0)
        % Received Signal (y = hs + n) Generation
        ReceivedSymbolSequence = H * SymbolSequence + NoiseSequence * sqrt(1 / EsN0(indx_EbN0)); % log2(M)x1 matrix
        
        % MLD Receiver
        [BitErrorCount_tmp, SignalErrorCount_tmp] = simulate_mld(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H);
        BitErrorCount(7,indx_EbN0) = BitErrorCount(7,indx_EbN0) + BitErrorCount_tmp;
        SignalErrorCount(7,indx_EbN0) = SignalErrorCount(7,indx_EbN0) + SignalErrorCount_tmp;
        
        % ZF
        [BitErrorCount_tmp, SignalErrorCount_tmp] = simulate_zf(ReceivedSymbolSequence, SignalSequence, SignalBinary, M, H);
        BitErrorCount(1,indx_EbN0) = BitErrorCount(1,indx_EbN0) + BitErrorCount_tmp;
        SignalErrorCount(1,indx_EbN0) = SignalErrorCount(1,indx_EbN0) + SignalErrorCount_tmp;
        
        % ZF - SIC
        [BitErrorCount_tmp, SignalErrorCount_tmp] = simulate_sic(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0(indx_EbN0), 'zf');
        BitErrorCount(2,indx_EbN0) = BitErrorCount(2,indx_EbN0) + BitErrorCount_tmp;
        SignalErrorCount(2,indx_EbN0) = SignalErrorCount(2,indx_EbN0) + SignalErrorCount_tmp;
        
        % ZF - OSIC
        [BitErrorCount_tmp, SignalErrorCount_tmp] = simulate_osic(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0(indx_EbN0), 'zf');
        BitErrorCount(5,indx_EbN0) = BitErrorCount(5,indx_EbN0) + BitErrorCount_tmp;
        SignalErrorCount(5,indx_EbN0) = SignalErrorCount(5,indx_EbN0) + SignalErrorCount_tmp;
        
        % MMSE - SIC
        [BitErrorCount_tmp, SignalErrorCount_tmp] = simulate_sic(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0(indx_EbN0), 'mmse');
        BitErrorCount(4,indx_EbN0) = BitErrorCount(4,indx_EbN0) + BitErrorCount_tmp;
        SignalErrorCount(4,indx_EbN0) = SignalErrorCount(4,indx_EbN0) + SignalErrorCount_tmp;
        
        % MMSE - OSIC
        [BitErrorCount_tmp, SignalErrorCount_tmp] = simulate_osic(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0(indx_EbN0), 'mmse');
        BitErrorCount(6,indx_EbN0) = BitErrorCount(6,indx_EbN0) + BitErrorCount_tmp;
        SignalErrorCount(6,indx_EbN0) = SignalErrorCount(6,indx_EbN0) + SignalErrorCount_tmp;
        
        % MMSE
        [BitErrorCount_tmp, SignalErrorCount_tmp] = simulate_mmse(ReceivedSymbolSequence, SignalSequence, SignalBinary, M, H, EsN0(indx_EbN0));
        BitErrorCount(3,indx_EbN0) = BitErrorCount(3,indx_EbN0) + BitErrorCount_tmp;
        SignalErrorCount(3,indx_EbN0) = SignalErrorCount(3,indx_EbN0) + SignalErrorCount_tmp;
    end
    if mod(iTotal-100, FivePercent)==0
        ElapsedTime = toc;
        EstimatedTime = (NumberIteration-iTotal)*ElapsedTime;
        disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/NumberIteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
    end
end

BER = BitErrorCount / (LengthBitSequence * NumberIteration);
SER = SignalErrorCount / (LengthSignalSequence * NumberIteration);

% Plot
BER_Title = sprintf("BER for %d-QAM %dX%d MIMO", M, Nt, Nr);
SER_Title = sprintf("SER for %d-QAM %dX%d MIMO", M, Nt, Nr);
x_axis = "Es/No (dB)";
<<<<<<< HEAD
legend_order = ["ZF", "SIC-ZF", "MMSE", "SIC-MMSE", "OSIC-ZF", "OSIC-MMSE", "MLD"];
myplot(EsN0_dB, BER, BER_Title, x_axis, "BER", legend_order);
ylim([10^(-4) 1])
myplot(EsN0_dB, SER, SER_Title, x_axis, "SER", legend_order);
ylim([10^(-4) 1])

legend_order = ["SIC", "OSIC", "ZF"];
myplot(EsN0_dB, BER, BER_Title, x_axis, "BER", legend_order);
myplot(EsN0_dB, SER, SER_Title, x_axis, "SER", legend_order);