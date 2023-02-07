function [BitErrorCount, SignalErrorCount] = simulate_osic(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0, ReceiverType)
    Nt = size(H,2);
    Nr = size(H,1);
    NormalizationFactor = sqrt(2/3*(M-1) * Nt);
    persistent alphabet
    if isempty(alphabet)
        alphabet = qammod([0:M-1], M) / NormalizationFactor;
    end
    snr = EsN0 / (NormalizationFactor^2);
    DetectedSignalSequence = zeros(Nt,1);
    HasValue = false(Nt,1);
    
    for ii = 1:Nt
        if strcmp(ReceiverType, 'zf')
            w = NormalizationFactor * pinv(H); % pinv(H) = inv(H' * H) * H'
        else
            % TODO: not sure about size(H,2) part
            w = NormalizationFactor * inv(H' * H + size(H,2) / EsN0 * eye(size(H,2))) * H';
        end
        wH = abs(w*H).^2;
        
        % TODO: there is definitely an error here. need fix
        sinr = snr*diag(wH)./(snr*(sum(wH,2) - diag(wH))+sum(abs(w).^2,2));
        [val,idx] = max(sinr);
        DetectedSymbol = w(idx, :) * ReceivedSymbolSequence;
        DetectedSignal = qamdemod(DetectedSymbol, M);
        
        OriginalIndex = get_original_index(HasValue, idx);
        DetectedSignalSequence(OriginalIndex, 1) = DetectedSignal;
        HasValue(OriginalIndex) = true;
        
        RemodulatedSignal = alphabet(DetectedSignal+1);
        ReceivedSymbolSequence = ReceivedSymbolSequence - H(:,idx) * RemodulatedSignal;
        H(:,idx) = []; % remove column
%         if ii ~= Nr
%             H = H(:, 2:size(H,2));
%         end
    end
    DetectedBinary = de2bi(DetectedSignalSequence, log2(M), 'left-msb');
    BitErrorCount = sum(SignalBinary~=DetectedBinary, 'all');
    SignalErrorCount = sum(SignalSequence~=DetectedSignalSequence, 'all');
end

function OriginalIndex = get_original_index(HasValue, idx)
    OriginalIndex = 0;
    while idx
        OriginalIndex = OriginalIndex + 1;
        if ~HasValue(OriginalIndex)
            idx = idx - 1;
        end
    end
end
