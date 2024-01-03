%% Description
% Modulation Diversity, but also possible in situations where Nt!=Nr

function BER = simulate_modulation_diversity_improved(eta, Nt, Nr, iteration, EbN0_dB)
    dbstop if error
    %% DEBGUG
    if ~exist('eta', 'var')
        Nt = 2;
        Nr = 4;
        eta = 2;
        EbN0_dB = 0;
        iteration = 10^3;
    end
%     n = 2;
%     eta = 2;
%     EbN0_dB = 0:5;
%     iteration = 10^4;
    NONOISE = false;
    %% BEGIN
    % eta/2 bits per dimension
    % to represent k bits, we need 2^k levels
    M = 2^(eta/2); % Modulation Order for one dimension
    
    EbN0 = db2pow(EbN0_dB);
    % EsN0 = EbN0 * bits per signal
    EsN0 = EbN0 * eta;
    EsN0_dB = pow2db(EsN0);
    
%     NormalizationFactor = sqrt(2/3*(M-1) * n);
    NormalizationFactor = sqrt(2/3*(M^2-1)*Nt);
    
    switch Nt
        case 2
            lambda = (1+sqrt(5))/2; % lambda = (1-sqrt(5))/2;

            a = 1/sqrt(1+lambda^2);
            b = lambda*a;
            R = [a -b;
                b  a];
        case 3
            phi = 2*cos(6/7*pi);
            alpha = (1+phi)/(1+phi+phi^2);
            beta = phi*alpha;
            gamma = -phi/(1+phi)*alpha;

            R = [alpha, beta, -gamma;
                beta, gamma, -alpha;
                gamma, alpha, -beta];
        otherwise
            error('n has to be 2 or 3');
    end

    %% Code Word generation
    
    %% Timer 
    FivePercent = ceil(iteration/20);
    
    
    %% Candidate Generation
    Candidates = zeros(2*Nt, Nt, M^(2*Nt^2));

    for ii = 0:M^(2*Nt^2)-1
        tmp = ii;
        count = 0;
        while tmp~=0
            Candidates(floor(count/Nt)+1, mod(count, Nt)+1, ii+1) = mod(tmp, M);
            tmp = floor(tmp/M);
            count = count + 1;
        end
    end
    
    CandidateSymbol = pammod(Candidates, M) / NormalizationFactor;
    
    CandidateSymbol = pagetranspose(pagemtimes(R,'none', CandidateSymbol, 'transpose'));
    
    %% Candidate Gen, restart
    for page=1:M^(2*Nt^2)
        for ii=2:Nt
            CandidateSymbol(:, ii, page) = circshift(CandidateSymbol(:,ii, page), ii-1);
        end
    end
    
    % Resize Page
    CandidateSymbol = reshape(CandidateSymbol, 2*Nt^2, 1, []);
    
    BEC = zeros(length(EsN0), 1);
    EuclideanDistance = zeros(size(CandidateSymbol, 3), 1);
    
    %% Validating Normalization
    TotalSum = zeros(size(CandidateSymbol,1), size(CandidateSymbol, 2));
    
    for page=1:M^(2*Nt^2)
        TotalSum(:,:) = TotalSum(:,:) + (CandidateSymbol(:,:,page).^2);
    end
    
    %% Simulation
    for iTotal = 1:iteration
        if mod(iTotal-100, FivePercent)==0
            tic
        end
        SignalSequence = randi([0 M-1], 2*Nt, Nt);
        
        SymbolSequence = pammod(SignalSequence, M) / NormalizationFactor;
        RotatedSymbol = (R*SymbolSequence')';
        for ii=2:Nt
            RotatedSymbol(:, ii) = circshift(RotatedSymbol(:,ii), ii-1);
        end
        x2_r = RotatedSymbol(:);
        
        H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) ./ sqrt(2); % slow fading
        H_r = [real(H), -imag(H);
                  imag(H), real(H)];
              
        TransmittedSignal = kron(eye(Nt), H_r) * x2_r;
        
        CandidateHX = pagemtimes(kron(eye(Nt), H_r), CandidateSymbol);

        Noise = randn(2*Nr*Nt, 1) / sqrt(2); % Nt is multiplied because Nt dictates the number of time slots
        if NONOISE
            Noise = zeros(size(Noise));
        end
        
        for idx = 1:length(EsN0)
            ReceivedSymbol = TransmittedSignal + Noise / sqrt(EsN0(idx));
            
            for page=1:size(CandidateSymbol,3)
                EuclideanDistance(page) = sum(abs(ReceivedSymbol - CandidateHX(:,:,page)).^2);
            end
            [~, IndexOfMin] = min(EuclideanDistance);
            
%             DetectedBinary = reshape(de2bi(IndexOfMin-1, log2(M)*Nt*2, 'left-msb'),log2(M),[])';
            BEC(idx) = BEC(idx) + sum(SignalSequence~=Candidates(:,:,IndexOfMin), 'all');
        end

        if mod(iTotal-100, FivePercent)==0
            ElapsedTime = toc;
            EstimatedTime = (iteration-iTotal)*ElapsedTime;
            disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/iteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
        end
    end
    BitsPerIteration = 2 * Nt^2 * log2(M);
    BER = BEC / (iteration * BitsPerIteration); % 4 symols per iteration; log2(M) bit per signal
end