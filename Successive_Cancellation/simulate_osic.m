function [BitErrorCount, SignalErrorCount] = simulate_osic(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0, ReceiverType)
    % IMPORTANT TODO: relative order to absolute order
    % if first idx is 4, numbers equal to or bigger than 4 is increased in
    % absolute order
    Nt = size(H,2);
    Nr = size(H,1);
    NormalizationFactor = sqrt(2/3*(M-1) * Nt); % size(H,1) = Nt
    snr = EsN0 / (NormalizationFactor^2);
    DetectedSignalSequence = zeros(Nt,1);
    for ii = 1:Nt
        if strcmp(ReceiverType, 'zf')
            w = NormalizationFactor * pinv(H); % pinv(H) = inv(H' * H) * H'
        else
            % TODO: not sure what to do about Nt; should it be decreased???
            w = NormalizationFactor * inv(H' * H + Nt / EsN0 * eye(size(H,2))) * H';
        end
        wH = abs(w*H).^2;
        % TODO: optimize sum of w; w*conj(w) or abs(w).^2
        sinr = snr*diag(wH)./(snr*(sum(wH,2) - diag(wH))+sum(abs(w).^2,2))
        [val,idx] = max(sinr)
        DetectedSymbol = w(idx, :) * ReceivedSymbolSequence;
        DetectedSignal = qamdemod(DetectedSymbol, M);
        DetectedSignalSequence(idx, 1) = DetectedSignal;
        RemodulatedSignal = qammod(DetectedSignal, M);
        ReceivedSymbolSequence = ReceivedSymbolSequence - H(:,idx) * RemodulatedSignal / NormalizationFactor;
        H(:,idx) = []; % remove column
%         if ii ~= Nr
%             H = H(:, 2:size(H,2));
%         end
    end
    DetectedBinary = de2bi(DetectedSignalSequence, log2(M), 'left-msb');
    BitErrorCount = sum(SignalBinary~=DetectedBinary, 'all');
    SignalErrorCount = sum(SignalSequence~=DetectedSignalSequence, 'all');
end