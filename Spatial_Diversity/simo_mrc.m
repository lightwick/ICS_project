function [BitErrorCount, SignalErrorCount] = simo_mrc(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H)
    Nt = size(H,2);
    Nr = size(H,1);
    assert(Nt==1, 'Nt is not 1')
    NormalizationFactor = sqrt(2/3*(M-1) * Nt); % size(H,1) = Nt
    
    y = ReceivedSymbolSequence;
    z = H'*y;
    
    DetectedSignal = qamdemod(z/norm(H,'fro')*NormalizationFactor, M);
    DetectedBinary = de2bi(DetectedSignal, log2(M), 'left-msb');
    
    SignalErrorCount = (DetectedSignal~=SignalSequence);
    BitErrorCount = sum(SignalBinary~=DetectedBinary, 'all');
end