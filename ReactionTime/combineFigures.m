% port data over to a new figure
fig_keep = figure(1);
ax_keep = gca;
hold on
%
fig_move = figure(2); % data from this figure will move to the other figure

color_move = getColorFromList(1,3); % color of data moved to other figure
% color_move = 'k';
offset = [0]; % shifts the x_data this much

line_move = findobj(fig_move,'type','line');
scatter_move = findobj(fig_move,'type','scatter');
patch_move = findobj(fig_move,'type','patch');
errorbar_move = findobj(fig_move,'type','errorbar');

% for p = 1:numel(patch_move)
%     new_handle = copyobj(patch_move(p),findobj(fig_keep,'type','axes'));
%     new_handle.FaceColor = color_move;
%     new_handle.EdgeColor = color_move;
%     new_handle.XData = patch_move(p).XData + offset;
% end

for s = 1:numel(scatter_move)
    new_handle = copyobj(scatter_move(s),findobj(fig_keep,'type','axes'));
    new_handle.MarkerFaceColor = color_move;
    new_handle.MarkerEdgeColor = 'none';
    new_handle.XData = scatter_move(s).XData + offset;
    
end

for l = 1:numel(line_move)
    if(strcmpi(line_move(l).LineStyle,'-'))
        
    else
        new_handle = copyobj(line_move(l),findobj(fig_keep,'type','axes'));
        new_handle.Color = color_move;
        new_handle.XData = line_move(l).XData + offset;
    end
end

for e = 1:numel(errorbar_move)
    new_handle = copyobj(errorbar_move(e),findobj(fig_keep,'type','axes'));
    new_handle.Color = color_move;
    new_handle.XData = errorbar_move(e).XData + offset;
end

%%
formatForLee(gcf)
set(gca,'FontSize',16)
xlabel('Amplitude on each electrode (\muA)')
ylabel('RT (s)')
xlim([0,38])
ylim([0.1,0.35])
%%

ax = gca;
ax.XAxis.MinorTickValues = [25:25:150,250];
% ax.XAxis.TickLabels = };



%% recolor a plot
fig_idx = 1;
color_move = getColorFromList(1,4);
line_move = findobj(fig_idx,'type','line');
scatter_move = findobj(fig_idx,'type','scatter');
errorbar_move = findobj(fig_idx,'type','errorbar');
patch_move = findobj(fig_idx,'type','patch');

offset = -0;

for s = 1:numel(scatter_move)
    new_handle = scatter_move(s);
    new_handle.MarkerFaceColor = color_move;
    new_handle.MarkerEdgeColor = 'none';
    new_handle.XData = scatter_move(s).XData + offset;
end

for l = 1:numel(line_move)
    new_handle = line_move(l);
    if(numel(new_handle.YData) > 1 && new_handle.YData(1) == new_handle.YData(2))
        % do nothing ... 
    else
        new_handle.Color = color_move;
        new_handle.XData = line_move(l).XData + offset;
        new_handle.MarkerSize = 24;
    end
end

for e = 1:numel(errorbar_move)
    new_handle = errorbar_move(e);
    new_handle.Color = color_move;
    new_handle.XData = errorbar_move(e).XData + offset;
end

