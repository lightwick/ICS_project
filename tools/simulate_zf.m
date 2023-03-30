function [BER, SER] = simulate_zf(M, Nr, Nt, iteration, EbN0_dB)
    EbN0 = db2pow(EbN0_dB);
    EsN0 = EbN0 * log2(M);
    EsN0_dB = pow2db(EsN0);
    
    %% DEBUG
    NormalizationFactor = sqrt(2/3*(M-1)*Nt);
    
    %% Timer
    FivePercent = ceil(iteration/20);
    
    %% START
    BEC = zeros(length(EsN0), 1);
    SEC = zeros(length(EsN0), 1);
    
    for iTotal = 1:iteration
        if mod(iTotal-100, FivePercent)==0
            tic
        end
        % Bit Generation
        SignalSequence = randi([0 M-1], Nt, 1);
        SignalBinary = de2bi(SignalSequence, log2(M), 'left-msb');
        SymbolSequence = qammod(SignalSequence, M) / NormalizationFactor;
        H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) ./ sqrt(2); % Receiver x Transmitter
        NoiseSequence = (randn(Nr, 1) + 1j * randn(Nr, 1)) / sqrt(2); % Noise (n) Generation
        
        for EsN0_idx = 1:length(EsN0)
            ReceivedSymbol = H * SymbolSequence + NoiseSequence * sqrt(1 / EsN0(EsN0_idx)); % log2(M)x1 matrix
            w_zf = pinv(H);
            DetectedSymbolSequence_ZF = w_zf * ReceivedSymbol;
            
            DetectedSequence = qamdemod(DetectedSymbolSequence_ZF*NormalizationFactor, M); % Detection
            DetectedBinary = de2bi(DetectedSequence, log2(M), 'left-msb');

            BEC(EsN0_idx) = BEC(EsN0_idx) + sum(SignalBinary~=DetectedBinary, 'all');
            SEC(EsN0_idx) = SEC(EsN0_idx) + sum(SignalSequence~=DetectedSequence, 'all');
        end
        
        if mod(iTotal-100, FivePercent)==0
            ElapsedTime = toc;
            EstimatedTime = (iteration-iTotal)*ElapsedTime;
            disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/iteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
        end
    end
    
    SER = SEC/(iteration*2);
    BER = BEC/(iteration*2*log2(M));
end