function [BitErrorCount, SignalErrorCount] = simulate_sic(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0, ReceiverType)
    Nt = size(H,2);
    Nr = size(H,1);
    NormalizationFactor = sqrt(2/3*(M-1) * Nt); % size(H,1) = Nt
    
    DetectedSignalSequence = zeros(Nt,1);
    for ii = 1:Nt
        if strcmp(ReceiverType, 'zf')
            w = NormalizationFactor * pinv(H); % pinv(H) = inv(H' * H) * H'
        else
            % TODO: not sure what to do about Nt; should it be decreased???
            w = NormalizationFactor * inv(H' * H + Nt / EsN0 * eye(size(H,2))) * H';
        end
        w = w(1,:);
        DetectedSymbol = w * ReceivedSymbolSequence;
        DetectedSignal = qamdemod(DetectedSymbol, M);
        DetectedSignalSequence(ii, 1) = DetectedSignal;
        RemodulatedSignal = qammod(DetectedSignal, M);
%         DetectedBinary = de2bi(DetectedSignalSequence, log2(M), 'left-msb');
        ReceivedSymbolSequence = ReceivedSymbolSequence - H(:,1) * RemodulatedSignal / NormalizationFactor;
        H(:,1) = []; % remove first column
%         if ii ~= Nr
%             H = H(:, 2:size(H,2));
%         end
    end
    DetectedBinary = de2bi(DetectedSignalSequence, log2(M), 'left-msb');
    BitErrorCount = sum(SignalBinary~=DetectedBinary, 'all');
    SignalErrorCount = sum(SignalSequence~=DetectedSignalSequence, 'all');
end