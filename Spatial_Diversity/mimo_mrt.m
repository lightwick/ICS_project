function [BitErrorCount, SignalErrorCount] = mimo_mrt(SymbolSequence, SignalSequence, Noise, SignalBinary,  M, H)
    % TODO: Normalization
    Nt = size(H,2);
    Nr = size(H,1);
    
    % S = svd(A) returns the singular values of matrix A in descending order.
    [U,S,V] = svd(H);
    w =V(:,1);
    g = U(:,1);
    
%     NormalizationFactor = sqrt(2/3*(M-1) * norm(w, "fro")^2);
    NormalizationFactor = sqrt(2/3*(M-1));
    TransmitSymbol = w * SymbolSequence / NormalizationFactor;
    ReceivedSymbol = H * TransmitSymbol + Noise;

    Processing = g' * ReceivedSymbol;
    
    DetectedSignal = qamdemod(Processing/S(1,1)*NormalizationFactor, M);
    DetectedBinary = de2bi(DetectedSignal, log2(M), 'left-msb');
    
    SignalErrorCount = sum(DetectedSignal~=SignalSequence, 'all');
    BitErrorCount = sum(SignalBinary~=DetectedBinary, 'all');
end