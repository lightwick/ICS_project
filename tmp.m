close all
clear
clc
% Simulation
M = 16
    NumberOfSignals = 10^2;
    LengthBitSequence = NumberOfSignals*log2(M); % log2(M) bits per signal
    
    NumberIteration = 10^3;
    
    Es = 1;
    Normalization_Factor = sqrt(2/3*(M-1));
    
    EbN0_dB = -2 : 2 : 20;
    EbN0 = 10 .^ (EbN0_dB / 10);
    
    % Calulate EsN0_dB
    % method 1; not sure how much of an error this makes or its effect
    EsN0 = EbN0 * log2(M);
    EsN0_dB = 10*log10(EsN0);
    % method 2
    EsN0_dB = EbN0_dB + 10*log10(log2(M)); % derived from EsN0 = EbN0 * log2(M);
    
    ErrorCount_ZF = zeros(1, length(EbN0_dB));
    ErrorCount_MMSE = zeros(1, length(EbN0_dB));
    ErrorCount_MLD = zeros(1, length(EbN0_dB));
    
    alphabet_letter = [-(sqrt(M)-1):2:sqrt(M)-1] / Normalization_Factor;
    alphabet = [];
    for ii = 1:length(alphabet_letter)
        for jj = 1:length(alphabet_letter)
            alphabet = [alphabet, alphabet_letter(ii) + j*alphabet_letter(jj)];
        end
    end
    
    for iTotal = 1 : NumberIteration
        BitSequence = randi([0 1], 1, LengthBitSequence); % Bit Generation (BitSequence = rand(1, LengthBitSequence) > 0.5;)
        SymbolSequence = qammod(BitSequence.', M, 'InputType', 'bit', 'UnitAveragePower', 1).';
        %avgPower = mean(abs(SymbolSequence).^2)
        NoiseSequence = (randn(1, length(SymbolSequence)) + 1j * randn(1, length(SymbolSequence))) / sqrt(2); % Noise (n) Generation
        H = (randn(1, length(SymbolSequence)) + 1j * randn(1, length(SymbolSequence))) ./ sqrt(2); % Channel (h) Generation
        for indx_EbN0 = 1 : length(EbN0)
            ReceivedSymbolSequence = H .* SymbolSequence + NoiseSequence * sqrt(1 / EsN0(indx_EbN0)); % Received Signal (y = s + n) Generation
    
            % ZF Receiver
            DetectionSymbolSequence_ZF = ReceivedSymbolSequence ./ H; % Detection (Zero-Forcing: y / h)
    
            % MMSE Receiver
            rho = EsN0(indx_EbN0);
            %w_mmse = conj(H) ./ (H.*conj(H)+1/rho);
            w_mmse = rho*(H'*H)^(-1)*H';
            pause
            DetectionSymbolSequence_MMSE = ReceivedSymbolSequence .* w_mmse;
    
            % MLD Receiver; ZF MLD 차이점????
            arg = (ones(length(alphabet),1) * ReceivedSymbolSequence) - (alphabet.' * H);
            arg = arg .* conj(arg);
            [val,idx] = min(arg);
            DetectionSymbolSequence_MLD = alphabet(idx); % TODO: could possibly simplify it more
    
            % Symbol Sequence -> Bit Sequence
            DetectionBitSequence_ZF = qamdemod(DetectionSymbolSequence_ZF.', M, 'OutputType', 'bit', 'UnitAveragePower', 1)'; % Detection
            % TODO: implementation of DetectionBitSequence_MMSE
            DetectionBitSequence_MMSE = qamdemod(DetectionSymbolSequence_MMSE.', M, 'OutputType', 'bit', 'UnitAveragePower', 1)'; % tmp value;
            DetectionBitSequence_MLD = qamdemod(DetectionSymbolSequence_MLD.', M, 'OutputType', 'bit', 'UnitAveragePower', 1)';
    
            ErrorCount_ZF(1, indx_EbN0) = ErrorCount_ZF(1, indx_EbN0) + sum(DetectionBitSequence_ZF~=BitSequence);
            ErrorCount_MMSE(1, indx_EbN0) = ErrorCount_MMSE(1, indx_EbN0) + sum(DetectionBitSequence_MMSE~=BitSequence);
            ErrorCount_MLD(1, indx_EbN0) = ErrorCount_MLD(1, indx_EbN0) + sum(DetectionBitSequence_MLD~=BitSequence);
        end
    end
    
    BER_Simulation_ZF = ErrorCount_ZF / (LengthBitSequence * NumberIteration);
    BER_Simulation_MMSE = ErrorCount_MMSE / (LengthBitSequence * NumberIteration);
    BER_Simulation_MLD = ErrorCount_MLD / (LengthBitSequence * NumberIteration);
    
    BER_Theory = berfading(EbN0_dB, 'qam', M, 1); % not sure if 'dataenc' needs to be specified; I don't even know what it does
    
    % Plot
    figure()
    semilogy(EbN0_dB, BER_Theory, 'r--'); % bin
    hold on
    semilogy(EbN0_dB, BER_Simulation_ZF, 'bo'); % gray
    semilogy(EbN0_dB, BER_Simulation_MMSE, 'bx');
    semilogy(EbN0_dB, BER_Simulation_MLD, 'b^');
    
    
    axis([-2 20 10^-5 0.5])
    grid on
    legend('Theory (Rayleigh)', 'ZF (Rayleigh)', 'MMSE (Rayleigh)', 'MLD (Rayleigh)');
    xlabel('Eb/No [dB]');
    ylabel('BER');
    title('BER for QAM (M='+string(M)+')');