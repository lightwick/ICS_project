function plotErgodicCapcity(Nr, Nt, EsN0_dB, Marker, Color, MarkerSize)
    iTotal = 10^4;
    ErgodicCapacity = zeros(length(EsN0_dB), 1);
    EsN0 = db2pow(EsN0_dB);
    
    for idx=1:length(EsN0_dB)
        Capacity = zeros(iTotal, 1);
        for iteration=1:iTotal
            H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) / sqrt(2); % Receiver x Transmitter
            C = log2(det(eye(Nr)+EsN0(idx)/Nt*H*H'));
            Capacity(iteration) = C;
        end
        ErgodicCapacity(idx) = real(mean(Capacity));
    end
    plot(EsN0_dB, ErgodicCapacity, Marker, 'Color', Color, 'MarkerSize', MarkerSize);
end