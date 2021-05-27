function [data_struct] = getDecayRateWrapper(data_struct, field_name, input_data)

    % initiate variables for each field
    data_struct.([field_name,'_decay_rates']) = zeros(numel(data_struct.stim_chan),12);
    data_struct.([field_name,'_is_responsive']) = zeros(numel(data_struct.stim_chan),12);
    data_struct.([field_name,'_baseline_counts']) = zeros(numel(data_struct.stim_chan),12);
    data_struct.([field_name,'_distance_from_stim']) = zeros(numel(data_struct.stim_chan),1);
    data_struct.([field_name,'_response_amp']) = zeros(numel(data_struct.stim_chan),12);
    data_struct.([field_name,'_response_to_each_pulse']) = cell(numel(data_struct.stim_chan),1);
    data_struct.([field_name,'_chan_stim']) = []; 
    data_struct.([field_name,'_chan_rec']) = []; 
    data_struct.([field_name,'_monkey']) = []; % 1 = han, 0 = duncan

    % get decay rate and amp data for each unit
    for u = 1:numel(data_struct.(field_name))
        if(strcmpi(data_struct.(field_name){u}.monkey, 'Han')==1)
            map_data = loadMapFile('R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp');
        elseif(strcmpi(data_struct.(field_name){u}.monkey, 'Duncan')==1)
            map_data = loadMapFile('R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp');
        end

        % stim window,baseline window,bin size,min mean baseline firing,
        % map_data, only use times when stimulation is occuring
        % (relevant for intermittent stimulation)
        decay_data = getDecayRate(data_struct.(field_name){u},[0,3900],[-2000,-200],...
            input_data.bin_size,input_data.min_rate*input_data.bin_size/1000,map_data,input_data.is_intermittent,...
            input_data.response_amp_time,input_data.response_amp_num_pulses,input_data.response_amp_pulse_window); 

        data_struct.([field_name,'_decay_rates'])(u,:) = decay_data.param_list(:,2)';
        data_struct.([field_name,'_intercept'])(u,:) = decay_data.param_list(:,1)';
        data_struct.([field_name,'_is_responsive'])(u,:) = decay_data.is_responsive';
        data_struct.([field_name,'_is_responsive_nonstim'])(u,:) = decay_data.is_responsive_nonstim';
        data_struct.([field_name,'_baseline_counts'])(u,:) = decay_data.baseline_counts';
        data_struct.([field_name,'_distance_from_stim'])(u) = decay_data.distance_from_stim_chan;
        data_struct.([field_name,'_response_amp'])(u,:) = decay_data.response_amp';
        data_struct.([field_name,'_response_to_each_pulse']){u} = decay_data.response_to_each_pulse;
        data_struct.([field_name,'_chan_stim'])(u,1) = data_struct.(field_name){u}.CHAN_LIST{1};
        data_struct.([field_name,'_chan_rec'])(u,1) = data_struct.(field_name){u}.CHAN_REC;
        data_struct.([field_name,'_monkey'])(u,1) = strcmpi(data_struct.(field_name){u}.monkey,'Han');% 1 = han, 0 = duncan
        data_struct.([field_name,'_spike_lat'])(u,:) = decay_data.spike_lat;
    end 

end

