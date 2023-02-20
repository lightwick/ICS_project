function [BitErrorCount, SignalErrorCount] = miso_alamouti(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H)
    Nt = size(H,2);
    Nr = size(H,1);
    assert(Nt==2 && Nr==1, "Need H of size 1x2")
    
    Augmented_H = [H; conj(H(1,2)) -conj(H(1,1))];
    
    NormalizationFactor = sqrt(2/3*(M-1) * Nt);
    
    y = [ReceivedSymbolSequence(1,1); conj(ReceivedSymbolSequence(2,1))];
    z = Augmented_H' * y;
    
    FrobSquared = H*H'; % equivalent to norm(H,'fro')^2
    DetectedSignal = qamdemod(z/FrobSquared*NormalizationFactor, M);
    DetectedBinary = de2bi(DetectedSignal, log2(M), 'left-msb');
    
    SignalErrorCount = sum(DetectedSignal~=SignalSequence, 'all');
    BitErrorCount = sum(SignalBinary~=DetectedBinary, 'all');
end