function [BitErrorCount, SignalErrorCount] = simulate_sic(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, ReceiverType)
    Nt = size(H,2);
    Nr = size(H,1);
    NormalizationFactor = sqrt(2/3*(M-1) * Nt); % size(H,1) = Nt
    
    DetectedSignalSequence = zeros(length(ReceivedSymbolSequence),1);
    for ii = 1:size(H,1)
        if strcmp(ReceiverType, 'zf')
            w = NormalizationFactor * pinv(H); % pinv(H) = inv(H' * H) * H'
        else
            w = NormalizationFactor * inv(H' * H + Nt / EsN0 * eye(Nt)) * H';
        end
        w = w(1,:);
        DetectedSymbol = w * ReceivedSymbolSequence;
        DetectedSignal = qamdemod(DetectedSymbol, M);
        DetectedSignalSequence(ii, 1) = DetectedSignal;
        RemodulatedSignal = qammod(DetectedSignal, M);
%         DetectedBinary = de2bi(DetectedSignalSequence, log2(M), 'left-msb');
        ReceivedSymbolSequence = ReceivedSymbolSequence - H(:,1) * RemodulatedSignal / NormalizationFactor;
        if ii ~= Nr
            H = H(:, 2:Nt);
        end
    end
    DetectedBinary = de2bi(DetectedSignalSequence, log2(M), 'left-msb');
    BitErrorCount = sum(SignalBinary~=DetectedBinary, 'all');
    SignalErrorCount = sum(SignalSequence~=DetectedSignalSequence, 'all');
end