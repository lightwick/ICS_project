function [EntireCapacity, ErgodicCapacity, OutageCapacity] = getCapcity_CK(Nr, Nt, EsN0_dB, iTotal)
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
            
            p = 1;
            [~,S,~] = svd(H);
            % TODO: Examine below implementation
            % Not quite sure if this is the correct implementation to get the eigenvalue.
            % First. I implemented this so we cam get eigenvalues of non-square matrices, but WHAT DOES THAT EVEN MEAN?? A NON SQUARE MATRIX EIGENVALUE??
            % Second. not sure if there are negative eigenvalues.
            if min(size(S))==1
                EigenValues = S(1,1)^2;
            else
                EigenValues = diag(S).^2;
            end
            
            while true
                r = rank(H);
                PowerDistribution = Nt/(r-p+1)*(1+1/EsN0(EsN0_idx)*sum(1./EigenValues, 'all')) - Nt/EsN0(EsN0_idx)*(1./EigenValues);
                if PowerDistribution(r-p+1)>=0
                    break;
                else
                    p = p+1;
                    EigenValues(length(EigenValues))=[];
                end
            end
            % TODO: Revise whether it should be Nt or Nt-p; logical antenna vs physical antenna
            Capacity = sum(log2(1+EsN0(EsN0_idx)/Nt * PowerDistribution.*EigenValues), 'all');
            EntireCapacity(iteration, EsN0_idx) = Capacity;
        end
        %% Ergodic Capcity
%         ErgodicCapacity(idx) = real(mean(EntireCapacity));
        
        %% Outage Capcity
        [cdf,x] = ecdf(EntireCapacity(:, EsN0_idx));
        [~, Outage_idx] = min(abs(cdf-0.1));
        OutageCapacity(EsN0_idx) = x(Outage_idx);
    end
    %% Ergodic Capcity
    ErgodicCapacity = real(mean(EntireCapacity))';
end