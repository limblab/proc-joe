% black = 20180510/20180511
% red = 20180515/20180516
% green = 20180521
% blue = 20180528

% port data over to a new figure

fig_orig = figure(14);
fig_new = figure(10);
ax_new = gca;
color_new = 'r';
for c = 1:numel(fig_orig.Children.Children)
    xData = fig_orig.Children.Children(c).XData;
    yData = fig_orig.Children.Children(c).YData;
    marker = fig_orig.Children.Children(c).Marker;
    markerSize = fig_orig.Children.Children(c).MarkerSize;
    lw = fig_orig.Children.Children(c).LineWidth;
    ls = fig_orig.Children.Children(c).LineStyle;
    hold on
    plot(ax_new,xData,yData,'marker',marker,'markersize',markerSize,'linewidth',lw,'linestyle',ls,'color',color_new);
    hold on
end