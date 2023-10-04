function BER = simulate_n4_modulation_diversity(eta, n, iteration, EbN0_dB)
    %% Current Limitations
    % It only works for when eta=2.
    % Wrote comment on where to change if needed for cases where eta~=2
    dbstop if error
    %% DEBGUG
    n = 4;
    eta = 2;
    EbN0_dB = 0;
    iteration = 1;
    assert(n==4, 'only works when n==4');
    NONOISE = true;

    %% BEGIN
    % eta/2 bits per dimension
    % to represent k bits, we need 2^k levels
    M = 2^(eta/2); % Modulation Order for one dimension
    
    EbN0 = db2pow(EbN0_dB);
    % EsN0 = EbN0 * bits per signal
    EsN0 = EbN0 * eta;
    EsN0_dB = pow2db(EsN0);
    
%     NormalizationFactor = sqrt(2/3*(M-1) * n);
    NormalizationFactor = sqrt(2/3*(M^2-1)*n);

    lambda_2 = (1+sqrt(5))/2;
    lambda = 1/10*(sqrt(2)+1)*sqrt(50+10*sqrt(2));
    U = sqrt(lambda_2^2 + lambda + lambda_2^2 * lambda^2) / lambda_2;

    a = 1/(U*sqrt(1+lambda_2^2));
    b = lambda_2/(U * sqrt(1+lambda_2^2));
    c = lambda / (U * lambda_2);
    d = lambda / U;

    R1 = [a b;
                -b a];
    R2 = [c d;
                -d c];
    R = [R1 -R2;
        R2 R1];
    
    %% Timer 
    % FivePercent = ceil(iteration/20);
    
    %% Candidate Generation
    
    BEC = zeros(length(EsN0), 1);
    EuclideanDistance = zeros(M^(2*n^2), 1);
    Candidates = zeros(2*n, n);
    
    %% Simulation
    for iTotal = 1:iteration
        % if mod(iTotal-100, FivePercent)==0
        %     tic
        % end
        SignalSequence = randi([0 M-1], 2*n, n);
        if eta~=2
            SymbolSequence = pammod(SignalSequence, M) / NormalizationFactor;
        else
            SymbolSequence = (SignalSequence-0.5) / NormalizationFactor;
        end
        RotatedSymbol = (R*SymbolSequence')';
        for ii=2:n
            RotatedSymbol(:, ii) = circshift(RotatedSymbol(:,ii), ii-1);
        end
        x2_r = RotatedSymbol(:);
        
        H = (randn(n, n) + 1j * randn(n, n)) ./ sqrt(2); % slow fading
        H_r = [real(H), -imag(H);
                  imag(H), real(H)];
              
        TransmittedSignal = kron(eye(n), H_r) * x2_r;

        if NONOISE
            Noise = zeros(2*n^2,1);
        else
            Noise = randn(2*n^2, 1) / sqrt(2);
        end
        
        for idx = 1:length(EsN0)
            ReceivedSymbol = TransmittedSignal + Noise / sqrt(EsN0(idx));
            
            for ii = 0:M^(2*n^2)-1
                if mod(ii+1, 10^6)==0
                    disp('toc'+string(ii));
                    toc
                end
                tmp = ii;
                count = 0;
                while tmp~=0
                    Candidates(floor(count/n)+1, mod(count, n)+1) = mod(tmp, M);
                    tmp = floor(tmp/M);
                    count = count + 1;
                end

                if eta~=2
                    CandidateSymbol = pammod(Candidates, M) / NormalizationFactor;
                else
                    CandidateSymbol = (Candidates - 0.5) / NormalizationFactor;
                end
                CandidateSymbol = (R * CandidateSymbol.').';
                
                %% Candidate Gen, restart
                for i2=2:n
                    CandidateSymbol(:, i2) = circshift(CandidateSymbol(:,i2), i2-1);
                end
                % Resize Page
                CandidateSymbol = reshape(CandidateSymbol, 2*n^2, 1);
                CandidateHX = kron(eye(n),H_r) * CandidateSymbol;

                EuclideanDistance(ii+1) = sum(abs(ReceivedSymbol - CandidateHX).^2);
                if mod(ii, 10^6)==0
                    disp('tic'+string(ii));
                    tic
                end
            end
            
            [~, IndexOfMin] = min(EuclideanDistance);
            
%             DetectedBinary = reshape(de2bi(IndexOfMin-1, log2(M)*Nt*2, 'left-msb'),log2(M),[])';
            DetectedSymbol = reshape(de2bi(IndexOfMin-1, 2*n^2), [], 2*n)'; % if eta~=2, then this would need to be changed
            BEC(idx) = BEC(idx) + sum(SignalSequence~=DetectedSymbol, 'all');
        end

        % if mod(iTotal-100, FivePercent)==0
        %     ElapsedTime = toc;
        %     EstimatedTime = (iteration-iTotal)*ElapsedTime;
        %     disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/iteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
        % end
    end
    BER = BEC / (iteration* 2*n^2 * log2(M)); % 4 symols per iteration; log2(M) bit per signal
end