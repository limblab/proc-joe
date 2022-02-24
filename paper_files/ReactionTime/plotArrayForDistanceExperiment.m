%% setup data for plotting
% need to provide electrodeList (cell array of grouped electrodes),
% colorList (color for each group of electrodes), and alphaList (supposed to let you pair
% groups)
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    map_data = loadMapFile(mapFileName);

    % electrodeList is typically loaded in
    alphaList = [];
    for i = 1:numel(electrodeList)
        alphaList(i) = 0.7+0.3*mod(i,2);
    end
    
    colorList = {};
    for i = 1:numel(electrodeList)
        color_temp = getColorFromList(2,ceil(i/2)-1);
        colorList{i} = ((1-alphaList(i))*1+(alphaList(i)*color_temp));
        colorList{i} = min(colorList{i},1);
    end
        
%% plot array map for distance experiment

    f=figure();
    f.Name = 'Han_distanceExp_arrayMap';
    hold on
    axis square
    for chan_idx = 1:96 % for every channel, determine how many electrodes to plot
        % then partition square based on that
        num_plot = 0;
        color_plot = [];
        alpha_plot = [];
        for elec_list_idx = 1:numel(electrodeList)
            if(~isempty(find(electrodeList{elec_list_idx} == chan_idx)))
                num_plot = num_plot + 1;
                color_plot(end+1,:) = colorList{elec_list_idx};
                alpha_plot(end+1,1) = alphaList(elec_list_idx);
            end
        end
        
        if(num_plot~=0)
            % get where to plot
             map_idx = find(map_data.chan == chan_idx);
             pos = [map_data.row(map_idx),11-map_data.col(map_idx)];
             
             for i = 1:num_plot
                % plot patch based on color_plot and alpha_plot
                x_data = pos(1)+[(i-1)/num_plot,i/num_plot,i/num_plot,(i-1)/num_plot];
                y_data = pos(2)+[0,0,1,1];
                patch(x_data,y_data,color_plot(i,:),'FaceAlpha',alpha_plot(i),'LineStyle','none')
                 
             end
             
             % plot rectangle around each electrode that we plotted on to
             % make it clear
             rectangle('Position',[pos(1),pos(2),1,1],'LineWidth',1.5)
        end
    end
    
    plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5)
    axis off
