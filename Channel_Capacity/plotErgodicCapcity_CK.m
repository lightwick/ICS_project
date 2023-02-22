function plotErgodicCapcity_CK(Nr, Nt, EsN0_dB, Marker, Color, MarkerSize)
% TODO: dimension compatability
    iTotal = 10^4;
    ErgodicCapacity = zeros(length(EsN0_dB), 1);
    EsN0 = db2pow(EsN0_dB);
    
    for idx=1:length(EsN0_dB)
        Capacity = zeros(iTotal, 1);
        for iteration=1:iTotal
            H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) / sqrt(2); % Receiver x Transmitter
            
            p = 1;
            [~,S,~] = svd(H);
            % TODO: Examine below implementation
            % Not quite sure if this is the correct implementation to get the eigenvalue.
            % First. I implemented this so we cam get eigenvalues of non-square matrices, but WHAT DOES THAT EVEN MEAN?? A NON SQUARE MATRIX EIGENVALUE??
            % Second. not sure if there are negative eigenvalues.
            EigenValues = diag(S).^2;
            while true
                r = rank(H);
                PowerDistribution = Nt/(r-p+1)*(1+1/EsN0(idx)*sum(1./EigenValues, 'all')) - Nt/EsN0(idx)*(1./EigenValues);
                if PowerDistribution(r-p+1)>=0
                    break;
                else
                    p = p+1;
                    EigenValues(length(EigenValues))=[];
                end
            end
            % TODO: Revise whether it should be Nt or Nt-p; logical antenna vs physical antenna
            C = sum(log2(1+EsN0(idx)/Nt * PowerDistribution.*EigenValues), 'all');
            Capacity(iteration) = C;
        end
        ErgodicCapacity(idx) = real(mean(Capacity));
    end
    plot(EsN0_dB, ErgodicCapacity, Marker, 'Color', Color, 'MarkerSize', MarkerSize);
end