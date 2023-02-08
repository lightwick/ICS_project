function fig = myplot(x, y, title_in, xlabel_in, ylabel_in, legend_in)
    colors = ["#0000FF", "#FF0000", "#4DBEEE", "#D95319", "#77AC30", "#EDB120", "#7E2F8E"];
    
    fig = figure();
    semilogy(x, y(1,:), '.--', 'Color', colors(1),'MarkerSize', 15);
    hold on
    for ii=2:size(y,1)
        semilogy(x, y(ii,:), '.--','Color', colors(ii), 'MarkerSize', 15);
    end
    ylabel(ylabel_in);
    title(title_in);
    grid on
    legend(legend_in);
    xlabel(xlabel_in);
end    