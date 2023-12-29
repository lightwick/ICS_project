function BER = simulate_orthogonal_mld(M , iteration, EbN0_dB)
    %% DEBUG
    % EbN0_dB = 0;
    % M = 4;
    % iteration = 10^5;
    % 
    % EbN0 = db2pow(EbN0_dB);
    % EsN0 = EbN0 * log2(M);
    %% Parameters and Setting
    Nt = 4;
    Nr = 4;
    Np = 2;
    TimeFrame = 2;
    EbN0 = db2pow(EbN0_dB);
    EsN0 = EbN0 * log2(M);
    
    %% Timer
    FivePercent = ceil(iteration/20);
    
    %% START
    Candidate_tmp = get_candidates(M, Np) / sqrt(Np);

    Candidates = zeros(Nt, TimeFrame, M^(Np*TimeFrame));
    count = 1;
    for i1=1:size(Candidate_tmp, 2)
        for i2 = 1:size(Candidate_tmp, 2)
            tmp = [Candidate_tmp(:,i1) Candidate_tmp(:,i2)];
            Candidates(:,:,count) = shift(Nt,Np,tmp);
            count = count+1;
        end
    end
    
    BEC = zeros(length(EsN0), 1);
    
    for iTotal = 1:iteration
        if mod(iTotal-100, FivePercent)==0
            tic
        end
        % Bit Generation
        SignalSequence = randi([0 M-1], Np, TimeFrame);
        SignalBinary = de2bi(SignalSequence, log2(M), 'left-msb');
        % SignalBinary = de2bi(SignalSequence, log2(M))

        SymbolSequence = qammod(SignalSequence, M, 'UnitAveragePower', true);
        % s_eff = zeros(Nt, TimeFrame);
        % s_eff([1:Np], 1) = SymbolSequence(:, 1);
        % s_eff([Np+1:Nt], 2) = SymbolSequence(:, 2);
        s_eff = shift(Nt, Np, SymbolSequence);

        H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) ./ sqrt(2); % Receiver x Transmitter
        NoiseSequence = (randn(Nr, TimeFrame) + 1j * randn(Nr, TimeFrame)) / sqrt(2); % Noise (n) Generation
        NoiseSequence = zeros(size(NoiseSequence));

        Hs_Candidate = pagemtimes(H, Candidates);

        for EsN0_idx = 1:length(EsN0)
            ReceivedSymbol = H * s_eff / sqrt(Np) + NoiseSequence * sqrt(1 / EsN0(EsN0_idx));

            % EuclideanDistance = abs(ReceivedSymbol * ones(1,M^Nt) - H*Candidates).^2;
            EuclideanDistance = sum(abs(ReceivedSymbol-Hs_Candidate).^2, [1 2]);
            [~, idx] = min(EuclideanDistance);

            % DetectedBinary = reshape(de2bi(idx-1, log2(M)*Nt, 'left-msb'),log2(M),[])';
            % DetectedSequence = bi2de(DetectedBinary, 'left-msb');
            % Separate each digit into individual numbers
            charArray = dec2base(idx-1, M, Np*TimeFrame);
            numArray = zeros(numel(charArray), 1);
            for i = 1:numel(charArray)
                numArray(i) = str2num(charArray(i));
            end
            DetectedBinary = de2bi(numArray, log2(M), 'left-msb');

            ErrorCount = sum(SignalBinary~=DetectedBinary, 'all');
            BEC(EsN0_idx) = BEC(EsN0_idx) + ErrorCount;
            % SEC(EsN0_idx) = SEC(EsN0_idx) + sum(SignalSequence~=DetectedSequence, 'all');
        end
        
        if mod(iTotal-100, FivePercent)==0
            ElapsedTime = toc;
            EstimatedTime = (iteration-iTotal)*ElapsedTime;
            disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/iteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
        end
    end
    
    % SER = SEC/(iteration*Nt);
    BER = BEC/(iteration*Nt*log2(M));
end

function Candidates = get_candidates(M, Nt)
    AllNumbers = de2bi([0:M^Nt-1], Nt*log2(M), 'left-msb');
    Candidates = zeros(M^Nt, Nt);
    for ii = 1 : M^Nt
        for jj = 1 : Nt
            Candidates(ii,jj) = bi2de(AllNumbers(ii,log2(M)*(jj-1)+1:log2(M)*jj), 'left-msb');
        end
    end
    Candidates = qammod(Candidates',M, 'UnitAveragePower', true);
end

function eff = shift(Nt, Np, symbol)
    eff = zeros(Nt, size(symbol,2));
    eff([1:Np], 1) = symbol(:, 1);
    eff([Np+1:Nt], 2) = symbol(:, 2);
end