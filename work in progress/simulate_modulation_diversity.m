function [BER, SER] = simulate_modulation_diversity(eta, n, iteration, EbN0_dB)
    
    %% DEBGUG
    n = 2;
    eta = 2;
    EbN0_dB = 0:5;
    iteration = 1;
    
    %% BEGIN
    % eta/2 bits per dimension = 4 levels
    % to represent k bits, we need 2^k levels
    M = 2^(eta/2);
    
    EbN0 = db2pow(EbN0_dB);
    EsN0 = EbN0 * log2(M);
    EsN0_dB = pow2db(EsN0);
    
%     NormalizationFactor = sqrt(2/3*(M-1) * n);
    NormalizationFactor = sqrt(n);
    
    if n==2
%         TODO: doesn't matter? or does matter?
        lambda = (1+sqrt(5))/2;
%         lambda = (1-sqrt(5))/2;
        a = 1/sqrt(1+lambda^2);
        b = lambda*a;
        R = [a -b;
               b  a];
    elseif n==3
        phi = 2*cos(6/7*pi);
        alpha = (1+phi)/(1+phi+phi^2);
        beta = phi*alpha;
        gamma = -phi/(1+phi)*alpha;
        
        R = [alpha, beta, -gamma;
               beta, gamma, -alpha;
               gamma, alhpa, -beta]
    else
        error('n has to be 2 or 3');
    end
    %% Code Word generation
    
    %% Timer 
    FivePercent = ceil(iteration/20);
    
    
    %% Candidate Generation
    Candidates = zeros(2*n, n, M^(2*n^2));
%     Candidates = zeros(M^Nt, Nt*2);

    for ii = 0:M^(2*n^2)-1
        tmp = ii;
        count = 0;
        while tmp~=0
            Candidates(floor(count/n)+1, mod(count, n)+1, ii+1) = mod(tmp, M);
            tmp = floor(tmp/M);
            count = count + 1;
        end
    end
    
    CandidateSymbol = pagetranspose(pagemtimes(R,'none', pammod(Candidates, M), 'transpose')) / NormalizationFactor
    
    for page=1:M^(2*n^2)
    for ii=2:2*n
            RotatedSymbol(:, ii) = circshift(RotatedSymbol(:,ii), ii+1);
    end
        
%      for ii = 0:M^(Nt*2)-1
%          for jj = 1:Nt*2
%              Candidates(ii+1,jj) = mod(floor(ii/M^(2*Nt-jj)),M);
%          end
%      end
%      Candidates = qammod(Candidates', M) / NormalizationFactor;
%      
%      power = zeros(1,2);
%     for ii = 1:M^(Nt*2)
%         a = Candidates(1,ii);
%         b = Candidates(2,ii);
%         c = Candidates(3,ii);
%         d = Candidates(4,ii);
%         X = 1/sqrt(5)*[alpha*(a+b*theta), alpha*(c+d*theta); 1j*(alpha_hat*(c+d*theta_hat)), alpha_hat*(a+b*theta_hat)];
%         power = power + sum(abs(X).^2);
%         CandidateSymbol(:,ii) = [real(X(1,1)); imag(X(1,1)); real(X(2,1)); imag(X(2,1)); real(X(1,2)); imag(X(1,2)); real(X(2,2)); imag(X(2,2))];
%     end
%     
%     power = power/(M^(Nt*2)) % from this, we can see that the power is correctly normalized
    
    BEC = zeros(length(EsN0), 1);
    SEC = zeros(length(EsN0), 1);
    
    %% Simulation
    for iTotal = 1:iteration
        if mod(iTotal-100, FivePercent)==0
            tic
        end
        SignalSequence = randi([0 M-1], 2*n, n);
%         BitSequence = de2bi(SignalSequence, log2(M)*n, 'left-msb');
        
        SymbolSequence = pammod(SignalSequence, M) / NormalizationFactor;
        RotatedSymbol = (R*SymbolSequence')';
        for ii=2:2*n
            RotatedSymbol(:, ii) = circshift(RotatedSymbol(:,ii), ii+1);
        end
        x2_r = RotatedSymbol(:);
        
        H = (randn(n, n) + 1j * randn(n, n)) ./ sqrt(2); % slow fading
        H_r = [real(H), -imag(H);
                  imag(H), imag(H)];
              
        TransmittedSignal = kron(eye(n), H_r) * x2_r;
        
        CandidateHX = pagemtimes(kron(eye(n), H_r), )
%         Noise = (randn(n, n) + 1j * randn(n, n)) ./ sqrt(2);
        Noise = zeros(n,n);
        
        for idx = 1:length(EsN0)
            ReceivedSymbol = TransmittedSignal + Noise / sqrt(EsN0(idx));
            EuclideanDistance = abs(ReceivedSymbol * ones(1,M^(Nt*2)) - H*CandidateSymbol).^2;
            [~, mini] = min(sum(EuclideanDistance, 1));
            
            DetectedBinary = reshape(de2bi(mini-1, log2(M)*Nt*2, 'left-msb'),log2(M),[])';
            DetectedSignal = bi2de(DetectedBinary, 'left-msb');
%             DetectedSignal = bi2de(DetectedBinary, 'left-msb');
            BEC(idx) = BEC(idx) + sum(DetectedBinary~=BitSequence, 'all');
            SEC(idx) = SEC(idx) + sum(DetectedSignal~=SignalSequence);
        end
        if mod(iTotal-100, FivePercent)==0
            ElapsedTime = toc;
            EstimatedTime = (iteration-iTotal)*ElapsedTime;
            disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/iteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
        end
    end
    BER = BEC / (iteration*4 * log2(M)); % 4 symols per iteration; log2(M) bit per signal
    SER = SEC/ (iteration*4);
end