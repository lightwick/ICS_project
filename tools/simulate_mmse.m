function [BitErrorCount, SignalErrorCount] = simulate_mmse(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0)
    Nt = size(H,2);
    NormalizationFactor = sqrt(2/3*(M-1) * Nt);

    w_mmse = NormalizationFactor * inv(H' * H + Nt / EsN0 * eye(Nt)) * H';
    DetectedSymbolSequence_MMSE = w_mmse * ReceivedSymbolSequence;

    DetectedSignalSequence_MMSE = qamdemod(DetectedSymbolSequence_MMSE, M);
    DetectedBinary_MMSE = de2bi(DetectedSignalSequence_MMSE, log2(M), 'left-msb');

    BitErrorCount = sum(SignalBinary~=DetectedBinary_MMSE, 'all');
    SignalErrorCount = sum(SignalSequence~=DetectedSignalSequence_MMSE, 'all');
end