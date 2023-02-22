function plotCapacityCDF(Nr, Nt, EsN0)
    colors = ["#0000FF", "#FF0000", "#4DBEEE", "#D95319", "#77AC30", "#EDB120", "#7E2F8E"];

    iTotal = 10^4;
    Capacity = zeros(iTotal,1);
    
    for EsN0_idx = 1:length(EsN0)
        for iteration=1:iTotal
            H = (randn(Nr, Nt) + 1j * randn(Nr, Nt)) / sqrt(2); % Receiver x Transmitter
            C = log2(det(eye(Nr)+EsN0(EsN0_idx)/Nt*H*H'));
            Capacity(iteration) = C;
        end
        [cdf,x] = ecdf(Capacity);
        plot(x,cdf, "Color", colors(EsN0_idx));
        hold on
        
        ErgodicCapacity = real(mean(Capacity));
        ErgodicLine = xline(ErgodicCapacity,'-',{'Ergodic Capacity'});
        
        [~, idx] = min(abs(cdf-0.1));
        OutageCapacity_10 = x(idx); 
        plot([0 OutageCapacity_10], [cdf(idx) cdf(idx)], "k");
        plot([OutageCapacity_10 OutageCapacity_10], [0 cdf(idx)], "k");
        
        % using annotation function would be cleaner
        text(OutageCapacity_10, cdf(idx)+0.025, '\downarrow', 'HorizontalAlignment', 'center');
        text(OutageCapacity_10, cdf(idx)+0.07, '10% outage capacity', 'HorizontalAlignment', 'center');
%         text(OutageCapacity_10, cdf(idx), '\downarrow 10% outage capacity')
    end
    grid on
    xlabel("Rate (bps/Hz)");
    ylabel("CDF");
    title("CDF of Channel Rate");
end