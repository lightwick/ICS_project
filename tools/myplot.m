function fig = myplot(x, y, title_in, xlabel_in, ylabel_in, legend_in)
    colors = ["#0000FF", "#FF0000", "#4DBEEE", "#D95319", "#77AC30", "#EDB120", "#7E2F8E"];
    Marker = [".", "pentagram", "hexagram", "*", "square", "^", "diamond"]
    MarkerSize = 8 * ones(1, length(colors));
    MarkerSize(1) = 15;

    fig = figure();
    for ii=1:size(y,1)
        semilogy(x, y(ii,:), '--','Color', colors(ii), 'Marker', Marker(ii), 'MarkerSize', MarkerSize(ii));
        hold on
    end
    ylabel(ylabel_in);
    title(title_in);
    grid on
    legend(legend_in);
    xlabel(xlabel_in);
end    
