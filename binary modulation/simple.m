close all
clear
clc

NumberIteration = 10^4;
LengthBitSequence = 10^2;

EbN0_dB = -2 : 1 : 20;
EbN0 = 10 .^ (EbN0_dB / 10);

ErrorCount_ZF = zeros(1, length(EbN0_dB));
ErrorCount_MMSE = zeros(1, length(EbN0_dB));
ErrorCount_MLD = zeros(1, length(EbN0_dB));


for iTotal = 1 : NumberIteration
    BitSequence = randi([0 1], 1, LengthBitSequence); % Bit Generation (BitSequence = rand(1, LengthBitSequence) > 0.5;)
    % SymbolSequence = qammod(BitSequence, 2, 'bin'); % 1 and -1    
    SymbolSequence = 2 * BitSequence - 1; % Symbol (s) Generation; consisting of 1 and -1
    NoiseSequence = (randn(1, length(SymbolSequence)) + 1j * randn(1, length(SymbolSequence))) ./ sqrt(2); % Noise (n) Generation
    H = (randn(1, length(SymbolSequence)) + 1j * randn(1, length(SymbolSequence))) ./ sqrt(2); % Channel (h) Generation
    for indx_EbN0 = 1 : length(EbN0)
        ReceivedSymbolSequence_H = H .* SymbolSequence + NoiseSequence * sqrt(1 / EbN0(indx_EbN0)); % Received Signal (y = hs + n) Generation
        
        % ZF Receiver
        DetectionSymbolSequence_ZF = ReceivedSymbolSequence_H ./ H; % Detection (Zero-Forcing: y / h)
        
        % MMSE Receiver
        rho = EbN0(indx_EbN0);
        w_mmse = conj(H) ./ (H.*conj(H)+1/rho);
        DetectionSymbolSequence_MMSE = ReceivedSymbolSequence_H .* w_mmse;
        
        % MLD Receiver
        arg = ([1 1]' * ReceivedSymbolSequence_H) - ([-1 1]' * H);
        arg = arg .* conj(arg);
        DetectionSymbolSequence_MLD = (arg(1, :)-arg(2, :)) > 0; % TODO: could possibly simplify it more
        
        % Symbol Sequence -> Bit Sequence
        DetectionBitSequence_ZF = real(DetectionSymbolSequence_ZF)>0;
        DetectionBitSequence_MMSE = real(DetectionSymbolSequence_MMSE)>0; % not sure if this is right
        DetectionBitSequence_MLD = DetectionSymbolSequence_MLD;
        
        ErrorCount_ZF(1, indx_EbN0) = ErrorCount_ZF(1, indx_EbN0) + biterr(DetectionBitSequence_ZF, BitSequence);
        ErrorCount_MMSE(1, indx_EbN0) = ErrorCount_MMSE(1, indx_EbN0) + biterr(DetectionBitSequence_MMSE, BitSequence);
        ErrorCount_MLD(1, indx_EbN0) = ErrorCount_MLD(1, indx_EbN0) + biterr(DetectionBitSequence_MLD, BitSequence);
    end
end

scatterplot(SymbolSequence)

BER_Simulation_ZF = ErrorCount_ZF / (LengthBitSequence * NumberIteration);
BER_Simulation_MMSE = ErrorCount_MMSE / (LengthBitSequence * NumberIteration);
BER_Simulation_MLD = ErrorCount_MLD / (LengthBitSequence * NumberIteration);

BER_Theory2_H = berfading(EbN0_dB, 'psk', 2, 1);

% Plot
figure()
semilogy(EbN0_dB, BER_Theory2_H, 'k--');
hold on
%hold on don't know what it does
semilogy(EbN0_dB, BER_Simulation_ZF, 'bo');
semilogy(EbN0_dB, BER_Simulation_MMSE, 'rx');
semilogy(EbN0_dB, BER_Simulation_MLD, 'g^');
axis([-2 20 10^-5 0.5]) % axis([a b c d]) x-axis from a to b, y-axis from c to d
grid on
legend('Theory', 'ZF', 'MMSE', 'MLD');
xlabel('Eb/No [dB]');
ylabel('BER');
title('BER for Binary Modulation');