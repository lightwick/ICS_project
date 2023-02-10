close all
clear all
clc
dbstop if error
dbstop if warning

addpath('../tools/')

% Environment Varible
M = 16
Nt = 4
Nr = 4
NumberIteration = 10^4;

% Simulation
LengthBitSequence = Nt * log2(M); % log2(M) bits per signal
LengthSignalSequence = Nt;

% EbN0_dB = 0:5:25;
% EbN0 = db2pow(EbN0_dB);
% 
% EsN0 = EbN0 * log2(M);
% EsN0_dB = pow2db(EsN0);

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

BEC_tmp = zeros(size(BitErrorCount));
SEC_tmp = zeros(size(SignalErrorCount));
    
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
        BEC = zeros(size(BitErrorCount,1), 1);
        SEC = zeros(size(SignalErrorCount,1), 1);
        % Received Signal (y = hs + n) Generation
        ReceivedSymbolSequence = H * SymbolSequence + NoiseSequence * sqrt(1 / EsN0(indx_EbN0)); % log2(M)x1 matrix
        
        % ZF
        [BEC(1), SEC(1)] = simulate_zf(ReceivedSymbolSequence, SignalSequence, SignalBinary, M, H);
        
        % ZF - SIC
        [BEC(2), SEC(2)] = simulate_sic(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0(indx_EbN0), 'zf');
        
        % MMSE
        [BEC(3), SEC(3)] = simulate_mmse(ReceivedSymbolSequence, SignalSequence, SignalBinary, M, H, EsN0(indx_EbN0));
                
        % MMSE - SIC
        [BEC(4), SEC(4)] = simulate_sic(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0(indx_EbN0), 'mmse');
        
        % ZF - OSIC
        [BEC(5), SEC(5)] = simulate_osic(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0(indx_EbN0), 'zf');

        % MMSE - OSIC
        [BEC(6), SEC(6)] = simulate_osic(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0(indx_EbN0), 'mmse');
                
        % MLD Receiver
        [BEC(7), SEC(7)] = simulate_mld(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H);
        
        BEC_tmp(:, indx_EbN0) = BEC;
        SEC_tmp(:, indx_EbN0) = SEC;
    end
    BitErrorCount = BitErrorCount + BEC_tmp;
    SignalErrorCount = SignalErrorCount + SEC_tmp;
    
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

legend_order = ["ZF", "SIC-ZF", "MMSE", "SIC-MMSE", "OSIC-ZF", "OSIC-MMSE", "MLD"];
myplot(EsN0_dB, BER, BER_Title, x_axis, "BER", legend_order);
ylim([10^(-4) 1])
myplot(EsN0_dB, SER, SER_Title, x_axis, "SER", legend_order);
ylim([10^(-4) 1])