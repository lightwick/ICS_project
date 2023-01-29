close all
clear
clc

Es = 1;

% Output Control
BER = true

% Environment Varible
M = 4
Nt = 4
Nr = Nt
NumberIteration = 10^4;
NumberOfSignals = 1;

% Simulation
% TODO: LengthBitSequence and Nt, Nr independent하게 구성
LengthBitSequence = Nt * NumberOfSignals *log2(M); % log2(M) bits per signal
NormalizationFactor = sqrt(2/3*(M-1)*Nt);

EsN0_dB = 0:5:50;
EsN0 = db2pow(EsN0_dB);

EbN0 = EsN0 / log2(M);
EbN0_dB = pow2db(EbN0);

BitErrorCount_ZF = zeros(1, length(EbN0_dB));
SignalErrorCount_ZF = zeros(1, length(EbN0_dB));

BitErrorCount_MLD = zeros(1, length(EbN0_dB));
SignalErrorCount_MLD = zeros(1, length(EbN0_dB));

BitErrorCount_MMSE = zeros(1, length(EbN0_dB));
SignalErrorCount_MMSE = zeros(1, length(EbN0_dB));


AllNumbers = de2bi([0:M^Nt-1], Nt*log2(M), 'left-msb');
Candidates = zeros(M^Nt, Nt);
for ii = 1 : M^Nt
    for jj = 1 : Nt
        Candidates(ii,jj) = bi2de(AllNumbers(ii,log2(M)*(jj-1)+1:log2(M)*jj), 'left-msb');
    end
end
Candidates = qammod(Candidates',M) / NormalizationFactor;

alphabet = qammod([0:M-1], M) / NormalizationFactor;
AvgPowerPerAntenna = mean(abs(alphabet) .^ 2)

FivePercent = ceil(NumberIteration/20);
for iTotal = 1 : NumberIteration
    if mod(iTotal-1, FivePercent)==0
        tic
    end
    % Bit Generation
    BitSequence = randi([0 M-1], Nt, 1);
    BitBinary = de2bi(BitSequence, log2(M), 'left-msb');
    SymbolSequence = qammod(BitSequence, M) / NormalizationFactor;
    
    NoiseSequence = (randn(Nt, 1) + 1j * randn(Nt, 1)) / sqrt(2); % Noise (n) Generation
    H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) ./ sqrt(2); % Receiver x Transmitter
    for indx_EbN0 = 1 : length(EbN0)
        % Received Signal (y = hs + n) Generation
        ReceivedSymbolSequence = H * SymbolSequence + NoiseSequence * sqrt(1 / EsN0(indx_EbN0)); % log2(M)x1 matrix
        
        % MLD Receiver
        EuclideanDistance = abs(ReceivedSymbolSequence * ones(1,M^Nt) - H*Candidates).^2; % results in Nt x M^Nt, each column representing each candidate symbol combination
        [val, idx] = min(sum(EuclideanDistance, 1));
        DetectedBinary_MLD = reshape(de2bi(idx-1, log2(M)*Nt, 'left-msb'),log2(M),[])'; % MOST LIKELY WRONG
        DetectedSequence_MLD = bi2de(DetectedBinary_MLD, 'left-msb');
        
        BitErrorCount_MLD(indx_EbN0) = BitErrorCount_MLD(indx_EbN0) + sum(BitBinary~=DetectedBinary_MLD, 'all');
        SignalErrorCount_MLD(indx_EbN0) = SignalErrorCount_MLD(indx_EbN0) + sum(BitSequence~=DetectedSequence_MLD, 'all');
        
        % ZF Receiver
        w_zf= inv(H' * H) * H'; % Moore-Penrose inverse
        DetectedSymbolSequence_ZF = w_zf * ReceivedSymbolSequence; % Detection (Zero-Forcing: y / h)
        
        DetectedSignalSequence_ZF = qamdemod(DetectedSymbolSequence_ZF*NormalizationFactor, M); % Detection
        DetectedBinary_ZF = de2bi(DetectedSignalSequence_ZF, log2(M), 'left-msb');
        
        BitErrorCount_ZF(indx_EbN0) = BitErrorCount_ZF(indx_EbN0) + sum(BitBinary~=DetectedBinary_ZF, 'all');
        SignalErrorCount_ZF(indx_EbN0) = SignalErrorCount_ZF(indx_EbN0) + sum(BitSequence~=DetectedSignalSequence_ZF, 'all');
        
        % MMSE Receiver
        w_mmse = inv(H' * H + Nt / EsN0(indx_EbN0) * eye(Nr)) * H';
        DetectedSymbolSequence_MMSE = w_mmse * ReceivedSymbolSequence; % Detection (Zero-Forcing: y / h)
        
        DetectedSignalSequence_MMSE = qamdemod(DetectedSymbolSequence_MMSE*NormalizationFactor, M); % Detection
        DetectedBinary_ZF = de2bi(DetectedSignalSequence_MMSE, log2(M), 'left-msb');
        
        BitErrorCount_MMSE(indx_EbN0) = BitErrorCount_MMSE(indx_EbN0) + sum(BitBinary~=DetectedBinary_ZF, 'all');
        SignalErrorCount_MMSE(indx_EbN0) = SignalErrorCount_MMSE(indx_EbN0) + sum(BitSequence~=DetectedSignalSequence_MMSE, 'all');
    end
    if mod(iTotal-1, FivePercent)==0
        ElapsedTime = toc;
        EstimatedTime = (NumberIteration-iTotal)*ElapsedTime;
        disp(sprintf("%.2f%%, estimated wait time %d minutes %d seconds", iTotal/NumberIteration*100, floor(EstimatedTime/60), round(mod(EstimatedTime, 60))))
    end
end

SER_MLD = SignalErrorCount_MLD / (LengthBitSequence * NumberIteration);
BER_MLD = BitErrorCount_MLD / (LengthBitSequence * NumberIteration);

SER_ZF = SignalErrorCount_ZF / (LengthBitSequence * NumberIteration);
BER_ZF = BitErrorCount_ZF / (LengthBitSequence * NumberIteration);

SER_MMSE = SignalErrorCount_MMSE / (LengthBitSequence * NumberIteration);
BER_MMSE = BitErrorCount_MMSE / (LengthBitSequence * NumberIteration);

% Plot
figure()

if BER_Output == true
%     semilogy(EsN0_dB, BER_MLD, 'o--','Color', '#D95319');
    semilogy(EsN0_dB, BER_MLD, '.--','Color', '#D95319', 'MarkerSize', 15);
    hold on
    semilogy(EsN0_dB, BER_ZF, 'b.--', 'MarkerSize', 15);
    semilogy(EsN0_dB, BER_MMSE, '.--', 'Color', "#4DBEEE", 'MarkerSize', 15);
else
    semilogy(EsN0_dB, SER_MLD, '.--', 'Color', '#D95319', 'MarkerSize', 15); % 주황
    hold on
    semilogy(EsN0_dB, SER_ZF, 'b.--', 'MarkerSize', 15);
    semilogy(EsN0_dB, SER_MMSE, '.--', 'Color', "#4DBEEE", 'MarkerSize', 15); 
end

grid on
legend('Theory (Rayleigh)', 'ZF (Rayleigh)', 'MMSE (Rayleigh)', 'MLD (Rayleigh)');
xlabel('Es/No [dB]');
ylabel('BER');
title('BER for QAM (M='+string(M)+')');