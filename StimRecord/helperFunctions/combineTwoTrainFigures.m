figure_base = figure(1);
figure_move = figure(2);

% moves plots from figure_move to figure_base, changes figure_move
% linestyle to '--'


for axis_idx = 1:numel(figure_base.Children)
    
    line_move = findobj(figure_move.Children(axis_idx),'type','line');
    line_base = findobj(figure_base.Children(axis_idx),'type','line');
    
    % check to make sure these are data lines and not the vertical lines
    for line_idx = 1:numel(line_move)
        line_base_color = line_base(line_idx).Color*1.2;
        line_base_color(line_base_color > 1) = 1;
        line_base(line_idx).Color = line_base_color;
        
        if(numel(line_move(line_idx).XData) > 2) % change linestyle, move
            line_move(line_idx).LineStyle = '-';
            line_move(line_idx).Color = line_move(line_idx).Color*0.55;
            new_handle = copyobj(line_move(line_idx),figure_base.Children(axis_idx));
            uistack(new_handle,'bottom');
        end
    end
    
end