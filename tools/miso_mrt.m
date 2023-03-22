function [BitErrorCount, SignalErrorCount] = miso_mrt(SymbolSequence, SignalSequence, Noise, SignalBinary,  M, H)
    % TODO: Normalization
    Nt = size(H,2);
    Nr = size(H,1);
    assert(Nt==2 && Nr==1, "Need H of size 1x2")

    NormalizationFactor = sqrt(2/3*(M-1) * norm(H, "fro")^2);
    w = H';
    TransmitSymbol = w * SymbolSequence / NormalizationFactor;
    ReceivedSymbol = H*TransmitSymbol + Noise;

    FrobSquared = H*H'; % equivalent to norm(H,'fro')^2
    DetectedSignal = qamdemod(ReceivedSymbol/FrobSquared*NormalizationFactor, M);
    DetectedBinary = de2bi(DetectedSignal, log2(M), 'left-msb');
    
    SignalErrorCount = sum(DetectedSignal~=SignalSequence, 'all');
    BitErrorCount = sum(SignalBinary~=DetectedBinary, 'all');
end