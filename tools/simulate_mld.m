function [BER, SER] = simulate_mld(M, n, iteration, EbN0_dB)
    EbN0 = db2pow(EbN0_dB);
    EsN0 = EbN0 * log2(M);
    EsN0_dB = pow2db(EsN0);
    
    %% DEBUG
    Nt = n;
    Nr = n;
    NormalizationFactor = sqrt(2/3*(M-1)*Nt);
    
    %% Timer
    FivePercent = ceil(iteration/20);
    
    %% START
    persistent Candidates
    if isempty(Candidates)
         Candidates = get_candidates(M, Nt) / NormalizationFactor;
    end
    
    BEC = zeros(length(EsN0), 1);
    SEC = zeros(length(EsN0), 1);
    
    
    for iTotal = 1:iteration
        if mod(iTotal-100, FivePercent)==0
            tic
        end
        % Bit Generation
        SignalSequence = randi([0 M-1], Nt, 1);
        SignalBinary = de2bi(SignalSequence, log2(M), 'left-msb');
        SymbolSequence = qammod(SignalSequence, M) / NormalizationFactor;
        H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) ./ sqrt(2); % Receiver x Transmitter
        NoiseSequence = (randn(Nr, 1) + 1j * randn(Nr, 1)) / sqrt(2); % Noise (n) Generation
        
        for EsN0_idx = 1:length(EsN0)
            ReceivedSymbol = H * SymbolSequence + NoiseSequence * sqrt(1 / EsN0(EsN0_idx)); % log2(M)x1 matrix
            % results in Nt x M^Nt, each column representing each candidate symbol combination
            EuclideanDistance = abs(ReceivedSymbol * ones(1,M^Nt) - H*Candidates).^2;
            [~, idx] = min(sum(EuclideanDistance, 1));

            DetectedBinary = reshape(de2bi(idx-1, log2(M)*Nt, 'left-msb'),log2(M),[])';
            DetectedSequence = bi2de(DetectedBinary, 'left-msb');

            BEC(EsN0_idx) = BEC(EsN0_idx) + sum(SignalBinary~=DetectedBinary, 'all');
            SEC(EsN0_idx) = SEC(EsN0_idx) + sum(SignalSequence~=DetectedSequence, 'all');
        end
        
        if mod(iTotal-100, FivePercent)==0
            ElapsedTime = toc;
            EstimatedTime = (iteration-iTotal)*ElapsedTime;
            disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/iteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
        end
    end
    
    SER = SEC/(iteration*2);
    BER = BEC/(iteration*2*log2(M));
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