%% plots a pattern of electrodes

% load map file
    input_data.mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    map_data = loadMapFile(input_data.mapFileName);
    
    chan_list = pattern1;
    
    input_data.color = getColorFromList(1,4);
%% plot heat map
%     figure();
    hold on
    for elec = 1:numel(chan_list)
        map_idx = find(chan_list(elec) == map_data.chan);
        rectangle('Position',[map_data.row(map_idx),map_data.col(map_idx),1,1],'FaceColor',input_data.color, 'EdgeColor','none');
        hold on
    end
    plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5)

    f=gcf;
    set(gca,'visible','off')
    xlim([1,11])
    ylim([1,11])
    axis square
