close all
clear
clc

NumberIteration = 10^4;
LengthBitSequence = 10^2;

Eb = 1;
EbN0_dB = -3 : 1 : 10;
EbN0 = 10 .^ (EbN0_dB / 10);

ErrorCount = zeros(1, length(EbN0_dB));
ErrorCount_H = zeros(1, length(EbN0_dB));

for iTotal = 1 : NumberIteration
    BitSequence = randi([0 1], 1, LengthBitSequence); % Bit Generation (BitSequence = rand(1, LengthBitSequence) > 0.5;)
    SymbolSequence = 2 * BitSequence - 1; % Symbol (s) Generation
    NoiseSequence = (randn(1, length(SymbolSequence)) + 1j * randn(1, length(SymbolSequence))) ./ sqrt(2); % Noise (n) Generation
    H = (randn(1, length(SymbolSequence)) + 1j * randn(1, length(SymbolSequence))) ./ sqrt(2); % Channel (h) Generation
    for indx_EbN0 = 1 : length(EbN0)
        ReceivedSymbolSequence = SymbolSequence + NoiseSequence * sqrt(1 / EbN0(indx_EbN0)); % Received Signal (y = s + n) Generation
        ReceivedSymbolSequence_H = H .* SymbolSequence + NoiseSequence * sqrt(1 / EbN0(indx_EbN0)); % Received Signal (y = hs + n) Generation
        DetectionSymbolSequence = real(ReceivedSymbolSequence) > 0; % Detection
        DetectionSymbolSequence_H = real(ReceivedSymbolSequence_H ./ H) > 0; % Detection (Zero-Forcing: y / h)
        DetectionBitSequence = DetectionSymbolSequence;
        DetectionBitSequence_H = DetectionSymbolSequence_H;
        ErrorCount_Tmp = sum(DetectionBitSequence ~= BitSequence); % Error Count
        ErrorCount_Tmp_H = sum(DetectionBitSequence_H ~= BitSequence); % Error Count
        ErrorCount(1, indx_EbN0) = ErrorCount(1, indx_EbN0) + ErrorCount_Tmp;
        ErrorCount_H(1, indx_EbN0) = ErrorCount_H(1, indx_EbN0) + ErrorCount_Tmp_H;
    end
end

BER_Simulation = ErrorCount / (LengthBitSequence * NumberIteration);
BER_Theory1 = 0.5 * erfc(sqrt(EbN0));
BER_Theory2 = berawgn(EbN0_dB, 'psk', 2, 'nondiff');

BER_Simulation_H = ErrorCount_H / (LengthBitSequence * NumberIteration);
BER_Theory2_H = berfading(EbN0_dB, 'psk', 2, 1);

% Plot
figure()
semilogy(EbN0_dB, BER_Simulation, 'r--');
hold on
semilogy(EbN0_dB, BER_Theory1, 'bo');
semilogy(EbN0_dB, BER_Theory2, 'bx');
semilogy(EbN0_dB, BER_Simulation_H, 'k--');
semilogy(EbN0_dB, BER_Theory2_H, 'kx');
axis([-3 10 10^-5 0.5])
grid on
legend('Simulation (AWGN)', 'Theory 1 (AWGN)', 'Theory 2 (AWGN)', 'Simulation (Fading)', 'Theory (Fading)');
xlabel('Eb/No [dB]');
ylabel('BER');
title('BER for Binary Modulation');