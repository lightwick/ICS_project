function [BitErrorCount, SignalErrorCount] = simulate_sic(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0, ReceiverType)
    Nt = size(H,2);
    Nr = size(H,1);
    NormalizationFactor = sqrt(2/3*(M-1) * Nt);
    persistent alphabet;
    if isempty(alphabet)
        alphabet = qammod([0:M-1], M) / NormalizationFactor;
    end
    DetectedSignalSequence = zeros(Nt,1);
    for ii = 1:Nt
        if strcmp(ReceiverType, 'zf')
            w = NormalizationFactor * pinv(H); % pinv(H) = inv(H' * H) * H'
        else
            w = NormalizationFactor * inv(H' * H + size(H,2) / EsN0 * eye(size(H,2))) * H';
        end
        w = w(1,:);
        DetectedSymbol = w * ReceivedSymbolSequence;
        DetectedSignal = qamdemod(DetectedSymbol, M);
        DetectedSignalSequence(ii, 1) = DetectedSignal;
        
        %% Remove the effect of the regarded transmit antenna
        RemodulatedSignal = alphabet(DetectedSignal+1);
        ReceivedSymbolSequence = ReceivedSymbolSequence - H(:,1) * RemodulatedSignal;
        H(:,1) = []; % remove first column
    end
    DetectedBinary = de2bi(DetectedSignalSequence, log2(M), 'left-msb');
    BitErrorCount = sum(SignalBinary~=DetectedBinary, 'all');
    SignalErrorCount = sum(SignalSequence~=DetectedSignalSequence, 'all');
end
