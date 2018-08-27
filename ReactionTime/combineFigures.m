% port data over to a new figure
fig_keep = figure(3);
ax_keep = gca;
fig_move = figure(1); % data from this figure will move to the other figure

color_move = getColorFromList(1,3); % color of data moved to other figure
offset = [-2.5]; % shifts the x_data this much

line_move = findobj(fig_move,'type','line');
scatter_move = findobj(fig_move,'type','scatter');
patch_move = findobj(fig_move,'type','patch');

for p = 1:numel(patch_move)
    patch_move(p).FaceColor = color_move;
    patch_move(p).EdgeColor = color_move;
    patch_move(p).XData = patch_move(p).XData + offset;
    copyobj(patch_move(p),findobj(fig_keep,'type','axes'));
end

for s = 1:numel(scatter_move)
    scatter_move(s).MarkerFaceColor = color_move;
    scatter_move(s).XData = scatter_move(s).XData + offset;
    copyobj(scatter_move(s),findobj(fig_keep,'type','axes'));
end

for l = 1:numel(line_move)
    line_move(l).Color = color_move;
    line_move(l).XData = line_move(l).XData + offset;
    copyobj(line_move(l),findobj(fig_keep,'type','axes'));
end

%%
formatForLee(gcf)
set(gca,'FontSize',16)
xlabel('Train length (ms)')
ylabel('RT (s)')
xlim([35,415])
ylim([0.1,0.35])

ax = gca;
ax.XAxis.MinorTickValues = [50:50:500];
