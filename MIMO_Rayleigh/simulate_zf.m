function [BitErrorCount, SignalErrorCount] = simulate_zf(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0)
        Nt = size(H,1);
        NormalizationFactor = sqrt(2/3*(M-1) * size(H,1)); % size(H,1) = Nt
        w_zf = (1/NormalizationFactor) * pinv(H); % pinv(H) = inv(H' * H) * H'
        DetectedSymbolSequence_ZF = w_zf * ReceivedSymbolSequence; % Detection (Zero-Forcing: y / h)
        
        DetectedSignalSequence_ZF = qamdemod(DetectedSymbolSequence_ZF, M); % Detection
        DetectedBinary_ZF = de2bi(DetectedSignalSequence_ZF, log2(M), 'left-msb');
        
        BitErrorCount = sum(SignalBinary~=DetectedBinary_ZF, 'all');
        SignalErrorCount = sum(SignalSequence~=DetectedSignalSequence_ZF, 'all');
end