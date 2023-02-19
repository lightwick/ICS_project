function [BitErrorCount, SignalErrorCount] = mimo_alamouti(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H)
% TODO: Review Normalization
    Nt = size(H,2);
    Nr = size(H,1);
    assert(Nt==2 && Nr==2, 'H is not size 2x2')
    assert(length(SignalSequence), "Signal Sequence is not 2")

    tmp_H = conj(H(:,[2 1]));
    tmp_H(:,2) = -tmp_H(:,2);

    Augmented_H = [H; tmp_H];
    NormalizationFactor = sqrt(2/3*(M-1) * Nt * 2);
    
    y = [ReceivedSymbolSequence(1,1); conj(ReceivedSymbolSequence(2,1))];
    z = Augmented_H' * y;
    
    FrobSquared = H*H'; % equivalent to norm(H,'fro')^2
    DetectedSignal = qamdemod(z/FrobSquared*NormalizationFactor, M);
    DetectedBinary = de2bi(DetectedSignal, log2(M), 'left-msb');
    
    SignalErrorCount = sum(DetectedSignal~=SignalSequence, 'all');
    BitErrorCount = sum(SignalBinary~=DetectedBinary, 'all');
end