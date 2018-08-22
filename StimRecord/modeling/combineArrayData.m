function [ data ] = combineArrayData( input_data )


    data = [];
    data.distance = [];
    data.p_val = [];
    data.response = [];
    data.amplitude = [];
    data.peak_latency = [];
    
    cd(input_data.folderpath);
    folders = dir();

    for f = 1:numel(folders)
        if(numel(folders(f).name) > 3 && isempty(strfind(folders(f).name,'Figures'))) % remove .., ., etc.
            cd(folders(f).name);
            array_data_filename = dir('*array*');
            load([input_data.folderpath,folders(f).name,'\',array_data_filename(1).name]);

            % need to get map_data
            monkey_name_idx = 0;
            for m = 1:numel(input_data.monkey_names)
                if(~isempty(strfind(folders(f).name,input_data.monkey_names{m})))
                    monkey_name_idx = m;
                end
            end
            map_data = loadMapFile(input_data.mapFileName{monkey_name_idx});

            % for each unit, get data relevant to X and Y
            for arr = 1:numel(arrayData)
                if(isfield(arrayData{arr},'STIM_PARAM_LIST')) % newer way of processing data
                    % chan stim, amplitude
                    for stim_param_idx = 1:size(arrayData{arr}.STIM_PARAM_LIST,1)
                        chan_stim_idx = find(map_data.chan == arrayData{arr}.STIM_PARAM_LIST(stim_param_idx,1));
                        chan_stim_pos = [map_data.row(chan_stim_idx), map_data.col(chan_stim_idx)];
                        
                        bin_idx_poststim = [find(arrayData{arr}.bE{stim_param_idx} >= input_data.poststim_window(1),1,'first'),...
                                find(arrayData{arr}.bE{stim_param_idx} >= input_data.poststim_window(2),1,'first')];
                        bin_idx_prestim = [find(arrayData{arr}.bE{stim_param_idx} >= input_data.prestim_window(1),1,'first'),...
                                find(arrayData{arr}.bE{stim_param_idx} >= input_data.prestim_window(2),1,'first')];
                            
                        [~,p] = kstest2(arrayData{arr}.bC{stim_param_idx}(bin_idx_prestim(1):bin_idx_prestim(2)),arrayData{arr}.bC{stim_param_idx}(bin_idx_poststim(1):bin_idx_poststim(2)));
                        if(p < input_data.max_p_val)
                            data.distance(end+1,1) = 400*sqrt((chan_stim_pos(1)-arrayData{arr}.ROW).^2 + (chan_stim_pos(2)-arrayData{arr}.COL).^2);
                            data.amplitude(end+1,1) = arrayData{arr}.STIM_PARAM_LIST(stim_param_idx,2);
                            data.p_val(end+1,1) = p;
                            
                            [~,max_idx] = max(arrayData{arr}.bC{stim_param_idx}(bin_idx_poststim(1):bin_idx_poststim(2)+5));
                            bin_size = mode(diff(arrayData{arr}.bE{stim_param_idx}));
                            zero_bin_idx = find(arrayData{arr}.bE{stim_param_idx} >= 0,1,'first');

                            % find correct bins
                            num_poststim_bins = diff(bin_idx_poststim);
                            data.response(end+1,1) = sum(arrayData{arr}.bC{stim_param_idx}(bin_idx_poststim(1):bin_idx_poststim(2))) - ...
                                mean(arrayData{arr}.bC{stim_param_idx}(bin_idx_prestim(1):bin_idx_prestim(2)))*num_poststim_bins;
                        
                            if(p >= input_data.p_val)
                                data.peak_latency(end+1,1) = -100;
                            else % compute peak latency, inhibition, etc.
                                data.peak_latency(end+1,1) = getPeakLatency(arrayData{arr}.bC{stim_param_idx},bin_idx_prestim,bin_idx_poststim,...
                                    mode(diff(arrayData{arr}.stim_param_idx)),zero_bin_idx,input_data);
                            end
                        
                        end
                        
                    end
                else % old default way
                    for chan_idx = 1:size(arrayData{arr}.bC,1)
                        for wave_idx = 1:size(arrayData{arr}.bC,2)
                            chan_stim_idx = find(map_data.chan == arrayData{arr}.CHAN_LIST(chan_idx,1));
                            chan_stim_pos = [11-map_data.row(chan_stim_idx), map_data.col(chan_stim_idx)];

                            % find correct bins
                            bin_idx_poststim = [find(arrayData{arr}.bE{chan_idx,wave_idx} >= input_data.poststim_window(1),1,'first'),...
                                find(arrayData{arr}.bE{chan_idx,wave_idx} >= input_data.poststim_window(2),1,'first')];
                            bin_idx_prestim = [find(arrayData{arr}.bE{chan_idx,wave_idx} >= input_data.prestim_window(1),1,'first'),...
                                find(arrayData{arr}.bE{chan_idx,wave_idx} >= input_data.prestim_window(2),1,'first')];
                            
                            [~,p] = kstest2(arrayData{arr}.bC{chan_idx,wave_idx}(bin_idx_prestim(1):bin_idx_prestim(2)),arrayData{arr}.bC{chan_idx,wave_idx}(bin_idx_poststim(1):bin_idx_poststim(2)));
                        
                            if(p < input_data.max_p_val)
                                data.distance(end+1,1) = 400*sqrt((chan_stim_pos(1)-arrayData{arr}.ROW).^2 + (chan_stim_pos(2)-arrayData{arr}.COL).^2);
                                data.amplitude(end+1,1) = arrayData{arr}.STIM_PARAMETERS(wave_idx).amp1;
                                data.p_val(end+1,1) = p;
                                
                                [~,max_idx] = max(arrayData{arr}.bC{chan_idx,wave_idx}(bin_idx_poststim(1):bin_idx_poststim(2)+5));
                                bin_size = mode(diff(arrayData{arr}.bE{chan_idx,wave_idx}));
                                zero_bin_idx = find(arrayData{arr}.bE{chan_idx,wave_idx} >= 0,1,'first');
                                                                
                                num_poststim_bins = diff(bin_idx_poststim);
                                data.response(end+1,1) = sum(arrayData{arr}.bC{chan_idx,wave_idx}(bin_idx_poststim(1):bin_idx_poststim(2))) - ...
                                    mean(arrayData{arr}.bC{chan_idx,wave_idx}(bin_idx_prestim(1):bin_idx_prestim(2)))*num_poststim_bins;
                                if(p >= input_data.p_val)
                                    data.peak_latency(end+1,1) = -100;
                                else % compute peak latency, inhibition, etc.
                                    data.peak_latency(end+1,1) = getPeakLatency(arrayData{arr}.bC{chan_idx,wave_idx},bin_idx_prestim,...
                                        bin_idx_poststim,zero_bin_idx,mode(diff(arrayData{arr}.bE{chan_idx,wave_idx})),input_data);
                                end
                            
                            end
                            
                            
                        end
                    end
                end

            end
        end

        cd(input_data.folderpath);
    end

    
    %% compute num_responsive and total_cells at each provided distance
    dists = unique(data.distance);
    data.num_responsive = zeros(numel(dists),1);
    data.total_cells = zeros(numel(dists),1);
    data.dists = dists;
    data.dist_markers = [];
    
    for d_idx = 1:numel(dists)
        dist_mask = data.distance == dists(d_idx);
        data.total_cells(d_idx) = sum(dist_mask);
        data.num_responsive(d_idx) = sum(data.p_val(dist_mask) < input_data.p_val);
    end

    
end



function [peak_latency] = getPeakLatency(data,bin_idx_prestim,bin_idx_poststim,zero_bin_idx,bin_size,input_data)

    [~,max_idx] = max(data(bin_idx_poststim(1):bin_idx_poststim(2)+10));
    if(data(max_idx + bin_idx_poststim(1)-1) > 2*std(data(bin_idx_prestim(1):bin_idx_prestim(2))) + mean(data(bin_idx_prestim(1):bin_idx_prestim(2))))
        peak_latency = (max_idx + bin_idx_poststim(1) - zero_bin_idx - 0.5)*bin_size - input_data.wave_duration;
    else
        peak_latency = -10;
    end

end