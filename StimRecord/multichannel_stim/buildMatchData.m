function [ match_data ] = buildMatchData( array_data, params )
    % array_data contains spike trial times for 1 neuron for each condition,
    % number of trials for each condition, a list of channels that 
    % were stimulated for each condition, and the channel that the neuron
    % was recorded on

    % match_data outputs for each multichannel stim:
    %       the probability of a spike in a window for that condition
    %       independence prediction for that set of channels 
    %       list of channels that were stimulated on and
    %           the probability of spike for each channel
    %       position of channels on array
    %       waveform sent
    
    
    
    % deal with params
    if nargin == 1, params = []; end
    
    
    % DEFAULT PARAMETERS
    post_stim_window_size = 3; % ms
    baseline_window = [-15,-5]; % ms, based on stim onset
    monkey_name = 'none';
    bootstrap = 1;
    num_boot = 1000;
    
    post_stim_window_bounds = [1.75,25]; % ms, when finding peak, look within this window
    post_stim_window_step_size = 0.1; % in ms
    %%%%%%%%%%%%%%%%%%%%%%
    if ~isfield(params,'monkey_name'), error('no monkey name provided, can not assign map file.'); end

    if nargin > 1 && ~isempty(params)
        assignParams(who,params); % overwrite parameters
    else
        params = struct();
    end
    
    helper_params.post_stim_window_size = post_stim_window_size;
    helper_params.baseline_window = baseline_window;
    helper_params.post_stim_window_bounds = post_stim_window_bounds;
    helper_params.post_stim_window_step_size = post_stim_window_step_size;
    helper_params.monkey_name = monkey_name;
       
    
    %%%%%%
    % call helper function to get relevant data
    sample_array_data.CHAN_LIST = array_data.CHAN_LIST;
    sample_array_data.CHAN_REC = array_data.CHAN_REC;
    sample_array_data.numStims = array_data.numStims;

    sample_array_data.spikeTrialTimes = array_data.spikeTrialTimes;
    
    match_data = buildMatchDataHelper(sample_array_data,helper_params);
    
    if(bootstrap)
        disp('bootstrapping data');
        for n = 1:num_boot-1
            if(mod(n,100) == 0)
                disp(num2str(n));
            end
            % resample spikes
            spike_times = []; % will be a list of spike times across conditions
            spike_cond = []; % 2xn matrix, (i,j) in arrayData.spikeTrialTimes
            for i = 1:size(array_data.spikeTrialTimes,1)
                for j = 1:size(array_data.spikeTrialTimes,2)
                    spike_times = [spike_times,array_data.spikeTrialTimes{i,j}];
                    this_cond = [i;j];
                    spike_cond = [spike_cond, this_cond*ones(size(array_data.spikeTrialTimes{i,j}))];
                end
            end
            
            [~,sample_idx] = datasample(spike_times,numel(spike_times),'Replace',true);
            sample_spikes = spike_times(sample_idx);
            sample_cond = spike_cond(:,sample_idx);
            sample_array_data.spikeTrialTimes = cell(size(array_data.spikeTrialTimes));
            
            % use idx's to break up into conditions again
            for i = 1:size(array_data.spikeTrialTimes,1)
                for j = 1:size(array_data.spikeTrialTimes,2)
                    sample_mask = sample_cond(1,:) == i & sample_cond(2,:) == j;
                    sample_array_data.spikeTrialTimes{i,j} = sample_spikes(sample_mask);
                end
            end
            
            
            % run helper function
            match_data_temp = buildMatchDataHelper(sample_array_data,helper_params);
            
            % store data in a nice format
            match_data.together(:,:,end+1) = match_data_temp.together;
            match_data.individual(:,:,end+1) = match_data_temp.individual;
            match_data.independence(:,:,end+1) = match_data_temp.independence;
            match_data.wave(:,:,end+1) = match_data_temp.wave;
            
        end
    end

end





%%%%%%%%%%%%%%%%
% helper function to actually get data


function [match_data] = buildMatchDataHelper(array_data,helper_params)

    %%%%%%%%%%%%%%%%%%%%%
    % unpack helper params
    post_stim_window_size = helper_params.post_stim_window_size;
    baseline_window = helper_params.baseline_window;
    post_stim_window_bounds = helper_params.post_stim_window_bounds;
    post_stim_window_step_size = helper_params.post_stim_window_step_size;
    monkey_name = helper_params.monkey_name;
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % load correct map file
    switch monkey_name
        case 'Han'
            map_data = loadMapFile('R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp');
        case 'Duncan'
            map_data = loadMapFile('R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp');
        otherwise
            error('monkey name does not match list of monkey files. Please add monkey and map file to function');
    end
    
    %%%%%%%%%%%%%%%%%%%%
    % setup output variables
    match_data.together = [];
    match_data.individual = [];
    match_data.independence = [];
    match_data.wave = [];
    match_data.pos = [];
    match_data.chans = [];
    match_data.chan_rec = array_data.CHAN_REC;
    
    match_data.chan_rec_pos = [11-map_data.row(find(map_data.chan == array_data.CHAN_REC)),map_data.col(find(map_data.chan == array_data.CHAN_REC))];

    
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    % get baseline number of spikes across all conditions
    
    spike_times = []; % will be a list of spike times across conditions
    spike_cond = []; % 2xn matrix, (i,j) in arrayData.spikeTrialTimes
    for i = 1:size(array_data.spikeTrialTimes,1)
        for j = 1:size(array_data.spikeTrialTimes,2)
            spike_times = [spike_times,array_data.spikeTrialTimes{i,j}];
            this_cond = [i;j];
            spike_cond = [spike_cond, this_cond*ones(size(array_data.spikeTrialTimes{i,j}))];
        end
    end
    
    
    baseline_mask = spike_times > baseline_window(1)/1000 & spike_times < baseline_window(2)/1000;
    
    baseline_num_spikes = sum(baseline_mask)/sum(sum(array_data.numStims))*post_stim_window_size/diff(baseline_window);
    
    
    %%%%%%%%%%%%%%%%%%%
    % for each condition compute the prob of a spike (number of spikes
    % above baseline)
    prob_spike_each_cond = zeros(size(array_data.numStims));
    window_center_used = zeros(size(array_data.numStims));
    for i = 1:size(array_data.numStims,1)
        for j = 1:size(array_data.numStims,2)
            
            % find window pos that maximizes number of spikes
            window_centers = (post_stim_window_bounds(1):post_stim_window_step_size:post_stim_window_bounds(2))'; % make a column vector
            windows = window_centers + [-post_stim_window_size/2,post_stim_window_size/2];
            
            spikes_window_mask = array_data.spikeTrialTimes{i,j} < windows(:,2)/1000 & array_data.spikeTrialTimes{i,j} > windows(:,1)/1000; 
            
            prob_spikes_window = sum(spikes_window_mask,2)/array_data.numStims(i,j);
            [~,max_idx] = max(prob_spikes_window); 
            prob_spike_each_cond(i,j) = prob_spikes_window(max_idx) - baseline_num_spikes;
            window_center_used(i,j) = window_centers(max_idx);
        end
    end
    
    %%%%%%%%%%%%%%%%% 
    % match individual electrode stimulation and group responses
    for i = 1:size(array_data.numStims,1)
        for j = 1:size(array_data.numStims,2)
             % if not a single electrode, store relevant data
            if(numel(array_data.CHAN_LIST{i}) > 1)
                match_data_idx = numel(match_data.wave)+1;
                match_data.pos{match_data_idx,1} = []; % [row col] for each elec (num stim elecs x 2 matrix)
                % find indices of single electrodes used
                elec_idx = [];
                for chan_idx = 1:numel(array_data.CHAN_LIST{i})
                    for list_idx = 1:numel(array_data.CHAN_LIST)
                        if(numel(array_data.CHAN_LIST{list_idx}) == 1 && array_data.CHAN_LIST{i}(chan_idx) == array_data.CHAN_LIST{list_idx})
                            elec_idx(end+1) = list_idx;
                            map_idx = find(map_data.chan == array_data.CHAN_LIST{elec_idx(end)});
                            % get elec position
                            match_data.pos{match_data_idx}(end+1,:) = [11-map_data.row(map_idx),map_data.col(map_idx)]; % [row col] for each elec (n x 2 matrix)
                        end
                    end                        
                end
                
                
                % store other data
                match_data.wave(match_data_idx,1) = j;
                
                match_data.together(match_data_idx,1) = prob_spike_each_cond(i,j);
                match_data.individual(match_data_idx,:) = prob_spike_each_cond(elec_idx,j)';
                match_data.independence(match_data_idx,1) = getIndependencePrediction(match_data.individual(match_data_idx,:));
                match_data.chans(match_data_idx,:) = [array_data.CHAN_LIST{elec_idx}];
                

            end % end if
            
        end
    end


end