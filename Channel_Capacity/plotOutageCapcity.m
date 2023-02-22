function plotOutageCapcity(Nr, Nt, EsN0_dB, Marker, Color, MarkerSize)
    iTotal = 10^4;
    OutageCapacity = zeros(length(EsN0_dB), 1);
    EsN0 = db2pow(EsN0_dB);
    
    for EsN0_idx=1:length(EsN0_dB)
        Capacity = zeros(iTotal, 1);
        for iteration=1:iTotal
            H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) / sqrt(2); % Receiver x Transmitter
            C = log2(det(eye(Nr)+EsN0(EsN0_idx)/Nt*H*H'));
            Capacity(iteration) = C;
        end
        [cdf,x] = ecdf(Capacity);
        [~, idx] = min(abs(cdf-0.1));
        OutageCapacity(EsN0_idx) = x(idx); 
    end
    plot(EsN0_dB, OutageCapacity, Marker, 'Color', Color, 'MarkerSize', MarkerSize);
end