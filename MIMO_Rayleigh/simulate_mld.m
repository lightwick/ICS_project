function [BitErrorCount, SignalErrorCount] = simulate_mld(ReceivedSymbolSequence, SignalSequence, SignalBinary,  M, H, EsN0)
    Nt = size(H,1);
    NormalizationFactor = sqrt(2/3*(M-1)*Nt);
    persistent Candidates
    if isempty(Candidates)
         Candidates = get_candidates(M, Nt) / NormalizationFactor;
    end
    
    EuclideanDistance = abs(ReceivedSymbolSequence/sqrt(EsN0) * ones(1,M^Nt) - H*Candidates).^2; % results in Nt x M^Nt, each column representing each candidate symbol combination
    [val, idx] = min(sum(EuclideanDistance, 1));
    
    DetectedBinary_MLD = reshape(de2bi(idx-1, log2(M)*Nt, 'left-msb'),log2(M),[])';
    DetectedSequence_MLD = bi2de(DetectedBinary_MLD, 'left-msb');

    BitErrorCount = sum(SignalBinary~=DetectedBinary_MLD, 'all');
    SignalErrorCount = sum(SignalSequence~=DetectedSequence_MLD, 'all');
end

function Candidates = get_candidates(M, Nt)
    AllNumbers = de2bi([0:M^Nt-1], Nt*log2(M), 'left-msb');
    Candidates = zeros(M^Nt, Nt);
    for ii = 1 : M^Nt
        for jj = 1 : Nt
            Candidates(ii,jj) = bi2de(AllNumbers(ii,log2(M)*(jj-1)+1:log2(M)*jj), 'left-msb');
        end
    end
    Candidates = qammod(Candidates',M);
end