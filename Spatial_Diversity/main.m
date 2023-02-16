close all
clear all
dbstop if error
% dbstop if warning

addpath('../tools/')

% Environment Varible
M = 4

NumberIteration = 10^4;
SimulationNum = 3;

% Simulation
EbN0_dB = 0:3:15;
EbN0 = db2pow(EbN0_dB);

EsN0 = EbN0 * log2(M);
EsN0_dB = pow2db(EsN0);

% EsN0_dB = 0:5:25;
% EsN0 = db2pow(EsN0_dB);
% 
% EbN0 = EsN0 / log2(M);
% EbN0_dB = pow2db(EbN0);

% BitErrorCount_ZF = zeros(1, length(EsN0_dB));
% SignalErrorCount_ZF = zeros(1, length(EsN0_dB));
% BitErrorCount_MLD = zeros(1, length(EsN0_dB));
% SignalErrorCount_MLD = zeros(1, length(EsN0_dB));
% 
% BitErrorCount_MMSE = zeros(1, length(EsN0_dB));
% SignalErrorCount_MMSE = zeros(1, length(EsN0_dB));

NormalizationFactor = sqrt(2/3*(M-1));

BitErrorCount = zeros(SimulationNum, length(EsN0_dB));
SignalErrorCount = zeros(SimulationNum, length(EsN0_dB));

BEC_tmp = zeros(size(BitErrorCount));
SEC_tmp = zeros(size(SignalErrorCount));
    
FivePercent = ceil(NumberIteration/20);
    
for iTotal = 1 : NumberIteration
    if mod(iTotal-100, FivePercent)==0
        tic
    end
    % Bit Generation
    SignalSequence = randi([0 M-1], 2, 1);
    SignalBinary = de2bi(SignalSequence, log2(M), 'left-msb');
    SymbolSequence = qammod(SignalSequence, M);
%     SymbolSequence = qammod(SignalSequence, M) / NormalizationFactor;
    
    NoiseSequence = (randn(2, 1) + 1j * randn(2, 1)) / sqrt(2); % Noise (n) Generation
%     NoiseSequence = zeros(2,1);
    H = (randn(2, 2) + 1j * randn(2, 2)) ./ sqrt(2); % Receiver x Transmitter

    for indx_EbN0 = 1 : length(EsN0)
        %% SIMO (MRC)
        Nt =1;
        Nr = 2;
        NormalizedSymbol = SymbolSequence(1:Nt, 1) / (NormalizationFactor * sqrt(Nt));
        y = H(1:Nr, 1:Nt) * NormalizedSymbol + NoiseSequence(1:Nt, 1) * sqrt(1 / EsN0(indx_EbN0));
        [BEC_tmp(1, indx_EbN0), SEC_tmp(1, indx_EbN0)] = simo_mrc(y(1:Nr,1), SignalSequence(1:Nt,1), SignalBinary(1:Nt, :),  M, H(1:Nr, 1:Nt));
        
        %% MISO (Alamouti)
        Nt = 2;
        Nr = 1;
        NormalizedSymbol = SymbolSequence(1:Nt, 1)/ (NormalizationFactor * sqrt(Nt * 2));
        % Row prepresents antenna, Column represents time-slot
        STBC = [NormalizedSymbol.'; -conj(NormalizedSymbol(2,1)) conj(NormalizedSymbol(1,1))].';
        y_alamouti = (H(1:Nr, 1:Nt) * STBC).' + NoiseSequence(1:Nt, 1) * sqrt(1 / EsN0(indx_EbN0));
        
        [BEC_tmp(2, indx_EbN0), SEC_tmp(2, indx_EbN0)] = miso_alamouti(y_alamouti, SignalSequence(1:Nt, 1), SignalBinary(1:Nt, :),  M, H(1:Nr,1:Nt));
        
        %% ZF (compare)
        Nt =2;
        Nr = 2;
        NormalizedSymbol = SymbolSequence(1:Nt, 1) / (NormalizationFactor * sqrt(Nt));
        y = H(1:Nr, 1:Nt) * NormalizedSymbol + NoiseSequence(1:Nt, 1) * sqrt(1 / EsN0(indx_EbN0));
        [BEC_tmp(3, indx_EbN0), SEC_tmp(3, indx_EbN0)] = simulate_zf(y, SignalSequence(1:Nt,1), SignalBinary(1:Nt, :),  M, H(1:Nr, 1:Nt));
    end
    
    BitErrorCount = BitErrorCount + BEC_tmp;
    SignalErrorCount = SignalErrorCount + SEC_tmp;
    
    if mod(iTotal-100, FivePercent)==0
        ElapsedTime = toc;
        EstimatedTime = (NumberIteration-iTotal)*ElapsedTime;
        disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/NumberIteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
    end
end

BER = BitErrorCount / (NumberIteration*log2(M)*2);
SER = SignalErrorCount / (NumberIteration*2);

BER(1,:) = BER(1,:) * 2;
SER(1,:) = SER(1,:) * 2;

% Plot
BER_Title = sprintf("BER for %d-QAM %dX%d MIMO", M, Nt, Nr);
SER_Title = sprintf("SER for %d-QAM %dX%d MIMO", M, Nt, Nr);
x_axis = "Eb/No (dB)";

legend_order = ["MRC", "Alamouti", "ZF"];
myplot(EbN0_dB, BER, BER_Title, x_axis, "BER", legend_order);
ylim([10^(-4) 1])
myplot(EbN0_dB, SER, SER_Title, x_axis, "SER", legend_order);
ylim([10^(-4) 1])

% colors = ["#0000FF", "#FF0000", "#4DBEEE", "#D95319", "#77AC30", "#EDB120", "#7E2F8E"];
% a=SER./BER;
% figure()
% plot(EsN0_dB, a(1,:), '.--', 'Color', colors(1), 'MarkerSize', 15);
% hold on
% for ii=2:size(a,1)
% plot(EsN0_dB, a(ii,:), '.--', 'Color', colors(ii), 'MarkerSize', 15);
% end
% ylabel("SER/BER")
% title("16-QAM 4x4 SER/BER")
% xlabel(x_axis)
% legend(legend_order)
% grid on
% title("16-QAM 2x2 SER/BER")