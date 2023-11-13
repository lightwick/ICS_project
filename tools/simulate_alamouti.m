function [BER, SER] = simulate_alamouti(M, Nr, iteration, EbN0_dB)
    if nargin==0
        M = 32;
        iteration = 10^6*5;
        Nr = 4;
        EsN0_dB = 0:2:22;
        EsN0 = db2pow(EsN0_dB);
        EbN0 = EsN0 / log2(M);
        EbN0_dB = pow2db(EbN0);
    end
    %% BEGIN
    EbN0 = db2pow(EbN0_dB);
    EsN0 = EbN0 * log2(M);
    EsN0_dB = pow2db(EsN0);
    
    Nt = 2;
  
    %% Timer
    FivePercent = ceil(iteration/20);
    
    BEC = zeros(length(EsN0), 1);
    SEC = zeros(length(EsN0), 1);
    
    %% Simulation
    for iTotal = 1:iteration
        if mod(iTotal-100, FivePercent)==0
            tic
        end
        SignalSequence = randi([0 M-1], 2, 1);
        BitSequence = de2bi(SignalSequence, log2(M), 'left-msb');
        
        
        H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) ./ sqrt(2); % Receiver x Transmitter
        Noise = (randn(Nr * 2, 1) + 1j * randn(Nr * 2, 1)) ./ sqrt(2); % Receiver x Transmitter
        
        for idx = 1:length(EsN0)
            SymbolSequence = qammod(SignalSequence, M, 'UnitAveragePower', true) * sqrt(EsN0(idx) / Nt);
            % STBC = [SymbolSequence.'; -conj(SymbolSequence(2,1)) conj(SymbolSequence(1,1))].';
            STBC = [SymbolSequence [-conj(SymbolSequence(2)) conj(SymbolSequence(1))].'];
            Hs = reshape((H * STBC), [], 1);
        
            y = Hs + Noise;
            
            Augmented_H = [H;
                                           conj(H(:, 2)) -conj(H(:,1))];

            y([Nr+1 : 2*Nr], :) = conj(y([Nr+1 : 2*Nr], :));
            z = Augmented_H' * y;
            FrobSquared = norm(H,'fro')^2;

            DetectedSignal = qamdemod(z / FrobSquared / sqrt(EsN0(idx)/Nt), M, 'UnitAveragePower', true);
            DetectedBinary = de2bi(DetectedSignal, log2(M), 'left-msb');

            BEC(idx) = BEC(idx) + sum(DetectedBinary~=BitSequence, 'all');
            SEC(idx) = SEC(idx) + sum(DetectedSignal~=SignalSequence);
        end
        if mod(iTotal-100, FivePercent)==0
            ElapsedTime = toc;
            EstimatedTime = (iteration-iTotal)*ElapsedTime;
            disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/iteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
        end
    end
    BER = BEC / (iteration* 2 * log2(M));
    SER = SEC / (iteration * 2);
end