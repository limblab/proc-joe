% black = 20180510/20180511
% red = 20180515/20180516
% green = 20180521
% blue = 20180528/20180531

% port data over to a new figure

fig_keep = figure(2); % data will be moved to this figure
fig_move = figure(3); % data from this figure will move to the other figure
color_move = getColorFromList(1,1); % color of data moved to other figure
offset = 0.1; % shifts the x_data this much

line_move = findobj(fig_move,'type','line');
scatter_move = findobj(fig_move,'type','scatter');

for l = 1:numel(line_move)
    line_move(l).Color = color_move;
    line_move(l).XData = line_move(l).XData + offset
    copyobj(line_move(l),findobj(fig_keep,'type','axes'));
end

for s = 1:numel(scatter_move)
    scatter_move(s).MarkerFaceColor = color_move;
    scatter_move(s).XData = scatter_move(s).XData + offset;
    copyobj(scatter_move(s),findobj(fig_keep,'type','axes'));
end