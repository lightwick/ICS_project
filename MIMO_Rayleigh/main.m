close all
clear
clc

Es = 1;

% Environment Varible
M = 4
Nt = 2
Nr = Nt
NumberIteration = 10^3;
SignalPerAntenna = 1;
NumberOfSignals = SignalPerAntenna;

% Simulation
% TODO: LengthBitSequence and Nt, Nr independent하게 구성
LengthBitSequence = Nt * SignalPerAntenna *log2(M); % log2(M) bits per signal
NormalizationFactor = sqrt(2/3*(M-1)) * Nt;

EsN0_dB = 0:5:25;
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
Candidates = Candidates';

alphabet = qamdemod(qammod([0:M-1], M),M);

for iTotal = 1 : NumberIteration
%     tic
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
        DetectedBinary_MLD = reshape(de2bi(idx-1, log2(M)*Nt, 'left-msb')',[],log2(M))' % MOST LIKELY WRONG, but still works on Nt=2
        DetectedBitSequence_MLD = bi2de(DetectedBinary_MLD, 'left-msb');
        
        BitErrorCount_MLD(indx_EbN0) = BitErrorCount_MLD(indx_EbN0) + sum(BitSequence~=DetectedBitSequence_MLD);
        SignalErrorCount_MLD(indx_EbN0) = SignalErrorCount_MLD(indx_EbN0) + sum(BitBinary~=DetectedBinary_MLD, 'all');
        
        % ZF Receiver
        w_zf= inv(H' * H) * H'; % Moore-Penrose inverse
        DetectedSymbolSequence_ZF = w_zf * ReceivedSymbolSequence; % Detection (Zero-Forcing: y / h)
        
        DetectedBitSequence_ZF = qamdemod(DetectedSymbolSequence_ZF*NormalizationFactor, M); % Detection
        DetectedBinary_ZF = de2bi(DetectedBitSequence_ZF, log2(M), 'left-msb');
        
        BitErrorCount_ZF(indx_EbN0) = BitErrorCount_ZF(indx_EbN0) + sum(BitSequence~=DetectedBitSequence_ZF);
        SignalErrorCount_ZF(indx_EbN0) = SignalErrorCount_ZF(indx_EbN0) + sum(BitBinary~=DetectedBinary_ZF, 'all');
%         ErrorCount_MMSE(1, indx_EbN0) = ErrorCount_MMSE(1, indx_EbN0) + sum(DetectionBitSequence_MMSE~=BitSequence);
%         ErrorCount_MLD(1, indx_EbN0) = ErrorCount_MLD(1, indx_EbN0) + sum(DetectionBitSequence_MLD~=BitSequence);
    end
%     toc
    if mod(iTotal, 10000)==0
        disp(string(iTotal/NumberIteration*100) + '% progress')
    end
end

SER_MLD = SignalErrorCount_MLD / (LengthBitSequence * NumberIteration);
BER_MLD = BitErrorCount_MLD / (LengthBitSequence * NumberIteration);

SER_ZF = SignalErrorCount_ZF / (LengthBitSequence * NumberIteration);
BER_ZF = BitErrorCount_ZF / (LengthBitSequence * NumberIteration);

SER_MMSE = SignalErrorCount_MMSE / (LengthBitSequence * NumberIteration);
BER_MMSE = BitErrorCount_MMSE / (LengthBitSequence * NumberIteration);
% BER_Simulation_ZF = ErrorCount_ZF / (LengthBitSequence * NumberIteration);
% BER_Simulation_MMSE = ErrorCount_MMSE / (LengthBitSequence * NumberIteration);
% BER_Simulation_MLD = ErrorCount_MLD / (LengthBitSequence * NumberIteration);

if M==2
    BER_Theory = berfading(EbN0_dB, 'psk', 2, 1);
else
    BER_Theory = berfading(EbN0_dB, 'qam', M, 1); % not sure if 'dataenc' needs to be specified; I don't even know what it does
end

% Plot
figure()
semilogy(EsN0_dB, BER_Theory, 'r--');
hold on
semilogy(EsN0_dB, SER_MLD, 'o', 'Color', '#D95319'); % 주황
semilogy(EsN0_dB, BER_MLD, 'x','Color', '#D95319');

semilogy(EsN0_dB, SER_ZF, 'bo');
semilogy(EsN0_dB, BER_ZF, 'bx');

semilogy(EsN0_dB, SER_ZF, 'bo');
semilogy(EsN0_dB, BER_ZF, 'bx');
% semilogy(EsN0_dB, BER_Simulation_MMSE, 'bx');
% semilogy(EsN0_dB, BER_Simulation_MLD, 'b^');


axis([0 25 10^-3 0.5])
grid on
% legend('Theory (Rayleigh)', 'ZF (Rayleigh)', 'MMSE (Rayleigh)', 'MLD (Rayleigh)');
xlabel('Es/No [dB]');
ylabel('BER');
title('BER for QAM (M='+string(M)+')');