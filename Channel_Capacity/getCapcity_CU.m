function [EntireCapacity, ErgodicCapacity, OutageCapacity] = getCapcity_CU(Nr, Nt, EsN0_dB, iTotal)
    % TODO: dimension compatability revision
    
    % EntireCapacity is a 'Column Vector'
    EntireCapacity = zeros(iTotal, length(EsN0_dB));
    ErgodicCapacity = zeros(length(EsN0_dB), 1);
    OutageCapacity = zeros(length(EsN0_dB), 1);
    
    EsN0 = db2pow(EsN0_dB);
    
    for EsN0_idx=1:length(EsN0_dB)
        %% Get Capacity array
        for iteration=1:iTotal
            H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) / sqrt(2); % Receiver x Transmitter
            Capacity = log2(det(eye(Nr)+EsN0(EsN0_idx)/Nt*H*H'));
            EntireCapacity(iteration, EsN0_idx) = Capacity;
        end
        
        %% Outage Capcity
        [cdf,x] = ecdf(EntireCapacity(:, EsN0_idx));
        [~, Outage_idx] = min(abs(cdf-0.1));
        OutageCapacity(EsN0_idx) = x(Outage_idx);
        
    end
    %% Ergodic Capcity
    ErgodicCapacity = real(mean(EntireCapacity))';
end