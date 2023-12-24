close all
clear all
dbstop if error
% dbstop if warning

addpath('../tools/')

% Environment Varible
M = 32;

NumberIteration = 10^5;
SimulationNum = 5;

% Simulation
% EbN0_dB = 0:3:15;
% EbN0 = db2pow(EbN0_dB);

% EsN0 = EbN0 * log2(M);
% EsN0_dB = pow2db(EsN0);

% EsN0_dB = 0:3:15;
% EsN0 = db2pow(EsN0_dB);
% 
% EbN0 = EsN0 / log2(M);
% EbN0_dB = pow2db(EbN0);

% BitErrorCount_ZF = zeros(1, length(EsN0_dB));
% SignalErrorCount_ZF = zeros(1, length(EsN0_dB));
% BitErrorCount_MLD = zeros(1, length(EsN0_dB));
% SignalErrorCount_MLD = zeros(1, length(EsN0_dB));
% 33edd
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
    
    NoiseSequence = (randn(4, 1) + 1j * randn(4, 1)) / sqrt(2); % Noise (n) Generation
%     NoiseSequence = zeros(4,1);
    H = (randn(2, 2) + 1j * randn(2, 2)) ./ sqrt(2); % Receiver x Transmitter

    for indx_EbN0 = 1 : length(EsN0)
        %% SIMO (MRC)
        Nt =1;
        Nr = 2;
        NormalizedSymbol = SymbolSequence(1:Nt, 1) / (NormalizationFactor * sqrt(Nt));
        Noise = NoiseSequence(1:Nr, 1) * sqrt(1 / EsN0(indx_EbN0));

        y = H(1:Nr, 1:Nt) * NormalizedSymbol + Noise;
        [BEC_tmp(1, indx_EbN0), SEC_tmp(1, indx_EbN0)] = simo_mrc(y, SignalSequence(1:Nt,1), SignalBinary(1:Nt, :),  M, H(1:Nr, 1:Nt));
        
        %% MISO (Alamouti)
        Nt = 2;
        Nr = 1;
        NormalizedSymbol = SymbolSequence(1:Nt, 1)/ (NormalizationFactor * sqrt(Nt));
        % Row prepresents antenna, Column represents time-slot
        STBC = [NormalizedSymbol.'; -conj(NormalizedSymbol(2,1)) conj(NormalizedSymbol(1,1))].';
        y_alamouti = (H(1:Nr, 1:Nt) * STBC).' + NoiseSequence(1:Nt, 1) * sqrt(1 / EsN0(indx_EbN0));
        
        [BEC_tmp(2, indx_EbN0), SEC_tmp(2, indx_EbN0)] = miso_alamouti(y_alamouti, SignalSequence(1:Nt, 1), SignalBinary(1:Nt, :),  M, H(1:Nr,1:Nt));
        
        %% MISO (MRT)
        Nt = 2;
        Nr = 1;
        H_new = H(1:Nr,1:Nt);
        [BEC_tmp(3, indx_EbN0), SEC_tmp(3, indx_EbN0)] = miso_mrt(SymbolSequence(1), SignalSequence(1), Noise(1:Nr), SignalBinary(1, :), M, H_new);

        %% MIMO (Alamouti)
        Nt = 2;
        Nr = 2;
        NormalizedSymbol = SymbolSequence(1:Nt, 1)/ (NormalizationFactor * sqrt(Nt));
        % Row prepresents antenna, Column represents time-slot
        STBC = [NormalizedSymbol.'; -conj(NormalizedSymbol(2,1)) conj(NormalizedSymbol(1,1))].';
        Hs = reshape((H(1:Nr, 1:Nt) * STBC), [], 1);
        y_alamouti = Hs + NoiseSequence * sqrt(1 / EsN0(indx_EbN0));
        [BEC_tmp(4, indx_EbN0), SEC_tmp(4, indx_EbN0)] = mimo_alamouti(y_alamouti, SignalSequence(1:Nt, 1), SignalBinary(1:Nt, :),  M, H(1:Nr, 1:Nt));
        
        %% MIMO (MRT)
        Nt = 2;
        Nr = 2;
        H_new = H(1:Nr,1:Nt);
        Noise = NoiseSequence(1:Nr, 1) * sqrt(1 / EsN0(indx_EbN0));
        [BEC_tmp(5, indx_EbN0), SEC_tmp(5, indx_EbN0)] = mimo_mrt(SymbolSequence(1), SignalSequence(1), Noise(1:Nr), SignalBinary(1, :), M, H_new);
    end
    
    BitErrorCount = BitErrorCount + BEC_tmp;
    SignalErrorCount = SignalErrorCount + SEC_tmp;
    
    if mod(iTotal-100, FivePercent)==0
        ElapsedTime = toc;
        EstimatedTime = (NumberIteration-iTotal)*ElapsedTime;
        disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/NumberIteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
    end
end

BER = BitErrorCount / (NumberIteration*log2(M));
SER = SignalErrorCount / (NumberIteration);

BER(2,:) = BER(2,:) / 2;
SER(2,:) = SER(2,:) / 2

BER(4,:) = BER(4,:) / 2;
SER(4,:) = SER(4,:) / 2;

% Plot
BER_Title = sprintf("BER for %d-QAM", M);
SER_Title = sprintf("SER for %d-QAM", M);
x_axis = "Eb/No (dB)";

legend_order = ["1x2 SIMO", "2x1 MISO [CU]", "2x1 MISO [CK]", "2x2 MIMO [CU]", "2x2 MIMO [CK]"];
myplot(EsN0_dB, BER, BER_Title, x_axis, "BER", legend_order);
ylim([10^(-6) 1])
myplot(EsN0_dB, SER, SER_Title, x_axis, "SER", legend_order);
ylim([10^(-6) 1])

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