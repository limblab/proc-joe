function [output_data] = plotResponseAmpVsDistance(array_data,input_data)

    % makes plots (1 for each wave type) of the responsive metric vs.
    % distance from stim chan across all neurons and stimulated channels
    
    % load in map data
    
    ax_list = [];
    figure_handle = figure();
    figure_handle.Position = [65.8000 307.4000 1.4552e+03 420.0000];
    
    dist_data_all = [];
    response_data_all = [];
    wave_idx_all = [];
    for arr_idx = 1:numel(array_data)
        if(strcmpi(array_data{arr_idx}.monkey,'Duncan')==1)
            map_file_name = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';      
        elseif(strcmpi(array_data{arr_idx}.monkey,'Han')==1)
            map_file_name = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
        else
            error('no map file for this monkey found');
        end
        
        map_data = loadMapFile(map_file_name);
        % get distances to each stim chan
        dist_data = zeros(size(array_data{arr_idx}.binCounts,1),1);

        rec_map_idx = find(map_data.chan == array_data{arr_idx}.CHAN_REC);
        rec_row = map_data.row(rec_map_idx); 
        rec_col = map_data.col(rec_map_idx);

        for chan = 1:size(array_data{arr_idx}.binCounts,1)
            if(iscell(array_data{arr_idx}.CHAN_LIST))
                stim_map_idx = find(map_data.chan == array_data{arr_idx}.CHAN_LIST{chan});
            else
                stim_map_idx = find(map_data.chan == array_data{arr_idx}.CHAN_LIST(chan));
            end
            stim_row = map_data.row(stim_map_idx);
            stim_col = map_data.col(stim_map_idx);
            dist_data(chan) = 400*sqrt((rec_row-stim_row)^2 + (rec_col-stim_col)^2);
        end
        
        % get response data for each waveform
        for wave = 1:size(array_data{arr_idx}.binCounts,2)
            plot_mask = dist_data > 0 & ~isnan(array_data{arr_idx}.response_amp(:,wave));
            if(input_data.use_responsive_only)
                plot_mask = array_data{arr_idx}.is_responsive(:,wave) & dist_data > 0;
            end
            
            dist_data_all = [dist_data_all; dist_data(plot_mask==1)];
            response_data_all = [response_data_all; array_data{arr_idx}.response_amp(plot_mask==1,wave)];
            wave_idx_all = [wave_idx_all; wave + zeros(sum(plot_mask == 1),1)];
                       
        end
    end
    
    % plot response data for each waveform
    ax_list = [];
    counter = 1;
    for wave = 1:4%size(array_data{arr_idx}.binCounts,2)
        
        dist_data_use = dist_data_all(wave_idx_all == wave);
        response_data_use = response_data_all(wave_idx_all == wave);
        
        if(input_data.use_boxplots)
            % bin distance data, no subplots yo
            boxplot_params.linewidth = input_data.linewidth;
            boxplot_params.box_width = input_data.box_width;
            boxplot_params.whisker_width = boxplot_params.box_width-0.02;
            boxplot_params.outlier_marker = input_data.outlier_marker;
            boxplot_params.outlier_marker_size = input_data.outlier_marker_size;
            boxplot_params.use_log_x_scale = 0;
            
            boxplot_params.outlier_color = getColorFromList(1,input_data.color_idx(wave));
            boxplot_params.median_color = getColorFromList(1,input_data.color_idx(wave));
            boxplot_params.box_color = getColorFromList(1,input_data.color_idx(wave));
            boxplot_params.whisker_color = getColorFromList(1,input_data.color_idx(wave));
            
            distance_bin_edges = (400:input_data.distance_bin_size:5200);
            distance_bin_centers = distance_bin_edges(1:end-1) + input_data.distance_bin_size/2;
            [~,~,distance_bin_map] = histcounts(dist_data_use,distance_bin_edges);
            
            for distance_idx = 1:numel(distance_bin_centers)
                distance_responses = response_data_use(distance_bin_map == distance_idx);
                if(numel(distance_responses) > 0)
                    boxplot_wrapper(distance_bin_centers(distance_idx)+input_data.dist_offset(wave),...
                        distance_responses,boxplot_params);
                end
            end
            
        else
            ax_list(counter) = subplot(1,4,counter);
            hold on
            plot(dist_data_use,response_data_use,'.','color','k');%getColorFromList(1,input_data.color_idx(wave)))
        end
        

        [fits{wave},gofs{wave}] = fit(dist_data_use,response_data_use,'a*exp(-x/b)','StartPoint',[0.5,1000]);
        
        x_data_fit = [0:50:5000];
        y_data_fit = feval(fits{wave},x_data_fit);
        if(input_data.use_boxplots) % colored line, or no fit....
%             plot(x_data_fit,y_data_fit,'--','linewidth',2,'color',getColorFromList(1,input_data.color_idx(wave)));
        else % black line
            plot(x_data_fit,y_data_fit,'r--','linewidth',2);
        end
        
        formatForLee(gcf)
        set(gca,'fontsize',14);
        if(wave == 1)
            xlabel('Distance (\mum)');
            ylabel('Response amp (spks/stim)');
        end
        counter = counter + 1;
    end
    if(~isempty(ax_list))
        linkaxes(ax_list,'xy');
    end
    xlim([-100,4600])
    
    
    
    output_data.fits = fits;
    output_data.gof = gofs;
    output_data.figure_handle = figure_handle;
end




























