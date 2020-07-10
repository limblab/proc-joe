function [output_data] = getResponseAndDistanceData(array_data, input_data)
% input_data contains
    % is_model : boolean
    % sub_baseline: boolean, do we subtract baseline
    
    
% output data contains response amplitude for each neuron, distance data,
% and a bunch of meta data to facilitate data analysis

    % setup data fields
    response_amp = []; 
    latency = []; latency_amp = []; latency_diam = []; latency_num_stims =[];
    dist_from_stim = [];
    monkey = []; % first letter of monkey's name, model = 'M';
    amp = [];
    unit_id = [];
    
    % will be populated for experiment data
    chan_rec = [];
    chan_stim = [];
    
    % will be populated for model data
    cell_id = [];
    clone_num = [];
    diam = [];
    
    
    % get map files. Currently implemented for Han and Duncan
    map_fpath = {'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp','R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp'};
    map_monkey = {'Han','Duncan'};
    
    for i_fpath = 1:numel(map_fpath)
        map_data{i_fpath} = loadMapFile(map_fpath{i_fpath});
    end
    % for each unit, compute response amp for each amp and chan
    for i_unit = 1:numel(array_data)
        % get baseline val
        if(input_data.sub_baseline == 1 && ~input_data.is_model)
            % convert to spikes per whatever bin size is put in
            baseline_val = array_data{i_unit}.baseline_fr*diff(input_data.spike_window);            
        else % set to 0
            baseline_val = 0;
        end
        
        map_data_idx = find(strcmpi(map_monkey,array_data{i_unit}.monkey)); % only useful for experiment data
        
        % for each channel and amplitude, get response amp
        for i_chan = 1:size(array_data{i_unit}.spikeTrialTimes,1)
            
            if(~input_data.is_model) % get stim_chan location
                if(iscell(array_data{i_unit}.CHAN_LIST))
                    stim_chan_idx = find(map_data{map_data_idx}.chan == array_data{i_unit}.CHAN_LIST{i_chan});
                else
                    stim_chan_idx = find(map_data{map_data_idx}.chan == array_data{i_unit}.CHAN_LIST(i_chan));
                end
                
                stim_chan_loc = [11-map_data{map_data_idx}.row(stim_chan_idx),map_data{map_data_idx}.col(stim_chan_idx)]; % 11- because that's how array data stores ROW and COL...very annoying but whatever
            end
            
            for i_amp = 1:size(array_data{i_unit}.spikeTrialTimes,2)
                % store meta data
                monkey(end+1,1) = array_data{i_unit}.monkey(1);
                unit_id(end+1,1) = i_unit;
                amp(end+1,1) = array_data{i_unit}.STIM_PARAMETERS(i_amp).amp1;
                
                if(input_data.is_model)
                    % distance is norm of loc
                    dist_from_stim(end+1,1) = norm(array_data{i_unit}.loc,2);
                    cell_id(end+1,1) = array_data{i_unit}.cell_id;
                    clone_num(end+1,1) = array_data{i_unit}.clone_num;
                    diam(end+1,1) = array_data{i_unit}.diam;
                else % not model
                    % get distance using stim_chan_loc
                    dist_from_stim(end+1,1) = 400*norm([array_data{i_unit}.ROW - stim_chan_loc(1), array_data{i_unit}.COL - stim_chan_loc(2)]);
                    % other variables from array_data
                    chan_rec(end+1,1) = array_data{i_unit}.CHAN_REC;
                    chan_stim(end+1,1) = map_data{map_data_idx}.chan(stim_chan_idx); % convoluted way to get it....
                end
                
                % compute response amp and store
                % just count spikes in spike_window and divide by total
                % number of stims
                
                spike_mask = array_data{i_unit}.spikeTrialTimes{i_chan,i_amp} > input_data.spike_window(1) & ...
                    array_data{i_unit}.spikeTrialTimes{i_chan,i_amp} <= input_data.spike_window(2);
                
                if(input_data.is_model && sum(spike_mask) > 1)
                    spike_mask(2:end) = 0; % remove doublets
                end
                response_amp(end+1,1) = sum(spike_mask)/array_data{i_unit}.numStims(i_chan,i_amp) - baseline_val;
                
                if(sum(spike_mask) > 0)
                    latency(end+1:end+sum(spike_mask))= array_data{i_unit}.spikeTrialTimes{i_chan,i_amp}(spike_mask);
                    latency_amp(end+1:end+sum(spike_mask)) = i_amp;
                    latency_num_stims(end+1:end+sum(spike_mask)) = array_data{i_unit}.numStims(i_chan,i_amp);
                    if(input_data.is_model)
                        latency_diam(end+1:end+sum(spike_mask)) = array_data{i_unit}.diam;
                    end
                end
            end
        end
    end

    % remove dists of 0
    if(~input_data.is_model)
        dist_mask = dist_from_stim <=0;
        response_amp(dist_mask) = [];
        dist_from_stim(dist_mask) = [];
        monkey(dist_mask) = [];
        amp(dist_mask) = [];
        unit_id(dist_mask) = [];
        chan_rec(dist_mask) = [];
        chan_stim(dist_mask) = [];
    end
    
    
    
    % package outputs
    output_data.response_amp = response_amp;
    output_data.dist_from_stim = dist_from_stim;
    output_data.monkey = monkey; % first letter of monkey's name, model = 'M';
    output_data.amp = amp;
    output_data.unit_id = unit_id;
    output_data.chan_rec = chan_rec;
    output_data.chan_stim = chan_stim;
    output_data.cell_id = cell_id;
    output_data.clone_num = clone_num;
    output_data.diam = diam;
    output_data.is_model = input_data.is_model;
    
    output_data.latency = latency;
    output_data.latency_amp = latency_amp;
    output_data.latency_diam = latency_diam;
    output_data.latency_num_stims = latency_num_stims;
end