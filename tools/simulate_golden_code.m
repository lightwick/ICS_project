function [BER, SER] = simulate_golden_code(M, iteration, EbN0_dB)
    %% BEGIN
    EbN0 = db2pow(EbN0_dB);
    EsN0 = EbN0 * log2(M);
    EsN0_dB = pow2db(EsN0);
    
    Nr = 2;
    Nt = 2;
    
    NormalizationFactor = sqrt(2/3*(M-1) * Nt);
    
    %% Code Word generation
    theta = (1+sqrt(5))/2;
    theta_hat = (1-sqrt(5))/2;
    alpha = 1+ 1j*(1-theta);
    alpha_hat = 1+1j*(1-theta_hat);
    
    %% Timer
    FivePercent = ceil(iteration/20);
    
    
    %% Candidate Generation
    Candidates = zeros(M^Nt, Nt*2);
    
     for ii = 0:M^(Nt*2)-1
         for jj = 1:Nt*2
             Candidates(ii+1,jj) = mod(floor(ii/M^(2*Nt-jj)),M);
         end
     end
     Candidates = qammod(Candidates', M) / NormalizationFactor;
     
     power = zeros(1,2);
    for ii = 1:M^(Nt*2)
        a = Candidates(1,ii);
        b = Candidates(2,ii);
        c = Candidates(3,ii);
        d = Candidates(4,ii);
        X = 1/sqrt(5)*[alpha*(a+b*theta), alpha*(c+d*theta); 1j*(alpha_hat*(c+d*theta_hat)), alpha_hat*(a+b*theta_hat)];
        power = power + sum(abs(X).^2);
        CandidateSymbol(:,ii) = [real(X(1,1)); imag(X(1,1)); real(X(2,1)); imag(X(2,1)); real(X(1,2)); imag(X(1,2)); real(X(2,2)); imag(X(2,2))];
    end
    
    power = power/(M^(Nt*2)) % from this, we can see that the power is correctly normalized
    
    BEC = zeros(length(EsN0), 1);
    SEC = zeros(length(EsN0), 1);
    
    %% Simulation
    for iTotal = 1:iteration
        if mod(iTotal-100, FivePercent)==0
            tic
        end
        SignalSequence = randi([0 M-1], 4, 1);
        BitSequence = de2bi(SignalSequence, log2(M), 'left-msb');
        
        SymbolSequence = qammod(SignalSequence, M) / NormalizationFactor;
        
        H = (randn(2, 2) + 1j * randn(2, 2)) ./ sqrt(2); % Receiver x Transmitter
        H = [real(H(1,1)), -imag(H(1,1)), real(H(1,2)), -imag(H(1,2));
            imag(H(1,1)), real(H(1,1)), imag(H(1,2)), real(H(1,2));
            real(H(2,1)), -imag(H(2,1)), real(H(2,2)), -imag(H(2,2));
            imag(H(2,1)), real(H(2,1)), imag(H(2,2)), real(H(2,2))];
        
        H = kron(eye(2), H);
        
        a = SymbolSequence(1);
        b = SymbolSequence(2);
        c = SymbolSequence(3);
        d = SymbolSequence(4);
        
        X = 1/sqrt(5)*[alpha*(a+b*theta), alpha*(c+d*theta); 1j*(alpha_hat*(c+d*theta_hat)), alpha_hat*(a+b*theta_hat)];
        X = [real(X(1,1)); imag(X(1,1)); real(X(2,1)); imag(X(2,1)); real(X(1,2)); imag(X(1,2)); real(X(2,2)); imag(X(2,2))];
        
        Noise = randn(8,1)/sqrt(2);
%         Noise = zeros(8,1);
        
        for idx = 1:length(EsN0)
            ReceivedSymbol = H*X + Noise / sqrt(EsN0(idx));
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