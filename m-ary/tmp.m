close all
clear
clc

% 수정사항: SNR이 클 때의 ZF, MMSE 비교

% Simulation
M = 16
Nt = 1;
NumberOfSignals = 1;
LengthBitSequence = Nt * NumberOfSignals*log2(M); % log2(M) bits per signal

NumberIteration = 10^4;

Es = 20:5:60;
% Normalization_Factor = sqrt(2/3*(M-1)); % 보고서에 해당 내용 정리

EsN0_dB = -2:2:20;
EsN0 = db2pow(EsN0_dB);

EbN0 = EsN0 / log2(M);
EbN0_dB = pow2db(EbN0);

ErrorCount_ZF = zeros(1, length(EbN0_dB));
ErrorCount_MMSE = zeros(1, length(EbN0_dB));
ErrorCount_MLD = zeros(1, length(EbN0_dB));

alphabet = qammod([0:M-1], M, 'UnitAveragePower', true);

for iTotal = 1 : NumberIteration
%     tic
    
    BitSequence = randi([0 1], 1, LengthBitSequence); % Bit Generation (BitSequence = rand(1, LengthBitSequence) > 0.5;)
    SymbolSequence = qammod(BitSequence.', M, 'InputType', 'bit', 'UnitAveragePower', 1).';
    %avgPower = mean(abs(SymbolSequence).^2)
    NoiseSequence = (randn(1, length(SymbolSequence)) + 1j * randn(1, length(SymbolSequence))) / sqrt(2); % Noise (n) Generation
    H = (randn(1, length(SymbolSequence)) + 1j * randn(1, length(SymbolSequence))) ./ sqrt(2); % Channel (h) Generation
    for indx_EbN0 = 1 : length(EbN0)
        ReceivedSymbolSequence = H .* SymbolSequence + NoiseSequence * sqrt(1 / EsN0(indx_EbN0)); % Received Signal (y = s + n) Generation

        % ZF Receiver
        w_zf = H.^(-1);
        DetectionSymbolSequence_ZF = ReceivedSymbolSequence .* w_zf; % Detection (Zero-Forcing: y / h)

        % MMSE Receiver
        w_mmse = (H.*conj(H)+1/EsN0(indx_EbN0)).^(-1) .* conj(H);
        z = ReceivedSymbolSequence .* w_mmse;
%         arg = (ones(length(alphabet),1) * z) - (alphabet.' * H .* w_mmse);
%         arg = arg .* conj(arg);
%         [val,idx] = min(arg);
%         DetectionSymbolSequence_MMSE = alphabet(idx);
        DetectionSymbolSequence_MMSE = z;
        
        % z로 qamdemod 사용

        % MLD Receiver
        arg = (ones(length(alphabet),1) * ReceivedSymbolSequence) - (alphabet.' * H);
        arg = abs(arg).^2;
        [val,idx] = min(arg);

        % Symbol Sequence -> Bit Sequence
        DetectionBitSequence_ZF = qamdemod(DetectionSymbolSequence_ZF.', M, 'OutputType', 'bit', 'UnitAveragePower', 1)'; % Detection
        DetectionBitSequence_MMSE = qamdemod(DetectionSymbolSequence_MMSE.', M, 'OutputType', 'bit', 'UnitAveragePower', 1)'; % tmp value;
        DetectionBitSequence_MLD = reshape(de2bi(idx-1, log2(M), 'left-msb')', 1, []);

        ErrorCount_ZF(1, indx_EbN0) = ErrorCount_ZF(1, indx_EbN0) + sum(DetectionBitSequence_ZF~=BitSequence);
        ErrorCount_MMSE(1, indx_EbN0) = ErrorCount_MMSE(1, indx_EbN0) + sum(DetectionBitSequence_MMSE~=BitSequence);
        ErrorCount_MLD(1, indx_EbN0) = ErrorCount_MLD(1, indx_EbN0) + sum(DetectionBitSequence_MLD~=BitSequence);
    end
%     toc
    
end

BER_Simulation_ZF = ErrorCount_ZF / (LengthBitSequence * NumberIteration);
BER_Simulation_MMSE = ErrorCount_MMSE / (LengthBitSequence * NumberIteration);
BER_Simulation_MLD = ErrorCount_MLD / (LengthBitSequence * NumberIteration);

if M==2
    BER_Theory = berfading(EbN0_dB, 'psk', 2, 1);
else
    BER_Theory = berfading(EbN0_dB, 'qam', M, 1); % not sure if 'dataenc' needs to be specified; I don't even know what it does
end

% Plot
figure()
semilogy(EsN0_dB, BER_Theory, 'r--');
hold on
semilogy(EsN0_dB, BER_Simulation_ZF, 'bo');
semilogy(EsN0_dB, BER_Simulation_MMSE, 'bx');
semilogy(EsN0_dB, BER_Simulation_MLD, 'b^');


axis([-2 20 10^-3 0.5])
grid on
legend('Theory (Rayleigh)', 'ZF (Rayleigh)', 'MMSE (Rayleigh)', 'MLD (Rayleigh)');
xlabel('Es/No [dB]');
ylabel('BER');
title('BER for QAM (M='+string(M)+')');