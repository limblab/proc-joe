function [  ] = plotPatternDataHeatmap( pattern_data, input_data )

    map_filename = input_data.map_filename(8:end); map_data = loadMapFile(map_filename);
    num_colors = input_data.num_colors;
    base_color = input_data.base_color;
    colors = [(0:num_colors-1)'/num_colors,(0:num_colors-1)'/num_colors,(0:num_colors-1)'/num_colors].*base_color;

    for pattern_idx = 1:numel(pattern_data.pattern)
        figure();
        hold on
        
        for chan_idx = 1:numel(pattern_data.pattern{pattern_idx}.chans)
            map_idx = find(pattern_data.pattern{pattern_idx}.chans(chan_idx) == map_data.chan);
            
            color_idx = max(1,min(size(colors,1),ceil(pattern_data.pattern{pattern_idx}.stim_norm(chan_idx)*size(colors,1))));

            rectangle('Position',[map_data.row(map_idx),map_data.col(map_idx),1,1],'FaceColor',colors(color_idx,:), 'EdgeColor','none');
            hold on
            
        end
        plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5)

        f=gcf;
        set(gca,'visible','off')
        xlim([1,11])
        ylim([1,11])
        axis square

        b = colorbar;
        colormap(colors);
        set(gca,'Visible','off');
        b.FontSize = 14;
        b.Ticks = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
        b.TickDirection = 'out';
        b.Ticks = [-1,-0.5,0,0.5,1];
        b.Ticks = [0,0.25,0.5,0.75,1];
        b.TickLabels = {'0','','0.5','','1'};
        b.Label.String = 'Normalized stimulation amplitude';
        b.Label.FontSize = 16;

    end
end

