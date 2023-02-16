function [BitErrorCount, SignalErrorCount] = simulate_zf(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H)
    Nt = size(H,2);
    NormalizationFactor = sqrt(2/3*(M-1) * Nt);
    w_zf = pinv(H); % pinv(H) = inv(H' * H) * H'
    DetectedSymbolSequence_ZF = w_zf * ReceivedSymbolSequence;

    DetectedSignalSequence_ZF = qamdemod(DetectedSymbolSequence_ZF*NormalizationFactor, M);
    DetectedBinary_ZF = de2bi(DetectedSignalSequence_ZF, log2(M), 'left-msb');

    BitErrorCount = sum(SignalBinary~=DetectedBinary_ZF, 'all');
    SignalErrorCount = sum(SignalSequence~=DetectedSignalSequence_ZF, 'all');
end