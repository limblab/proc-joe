function [ output_data ] = getDecayRate( array_data, stim_window, baseline_window, bin_size, min_baseline_counts, map_data, is_intermittent_stim, response_amp_time,response_amp_num_pulses,spike_window)
% send in a single array_data and a bin size, this then fits with a double
% exponential and outputs the decay rate from that

    array_data = rebinArrayData(array_data,bin_size);
        
    % get index in binCounts from window
    stim_window_idx = [find(array_data.binEdges{1} >= stim_window(1),1,'first'),find(array_data.binEdges{1} >= stim_window(2),1,'first')];
    baseline_window_idx = [find(array_data.binEdges{1} >= baseline_window(1),1,'first'),find(array_data.binEdges{1} >= baseline_window(2),1,'first')];
    
    % for each condition, fit data
    fits = cell(size(array_data.binCounts));
    gof = cell(size(fits));
    param_list = zeros(numel(fits),2);
    is_responsive = zeros(numel(fits),1);
    is_responsive_nonstim = zeros(numel(fits),1);
    baseline_counts = zeros(numel(fits),1);
    response_to_each_pulse = cell(numel(fits),1);
    lat_of_spikes = cell(numel(fits),1);
    response_amp = zeros(numel(fits),1);
    
    base_spike_times = [];
    for condition = 1:numel(fits)
        x_data_baseline = array_data.binEdges{condition}(baseline_window_idx(1):baseline_window_idx(2))'/1000; % in seconds
        y_data_baseline = array_data.binCounts{condition}(baseline_window_idx(1):baseline_window_idx(2))';
        
        baseline_counts(condition) = mean(y_data_baseline);
        base_spike_times = [base_spike_times,...
            array_data.spikeTrialTimes{condition}(array_data.spikeTrialTimes{condition} >= baseline_window(1)/1000 & array_data.spikeTrialTimes{condition} <= baseline_window(2)/1000)];
    end
    
    boot_data = [];
    boot_data.n = 1000;
    boot_data.base_spike_times = base_spike_times;
    boot_data.bin_edges = array_data.binEdges{1}(baseline_window_idx(1)):diff(spike_window):array_data.binEdges{1}(baseline_window_idx(2));
    boot_data.bin_edges = boot_data.bin_edges/1000;
    boot_data.n_trials = 8*12;
    baseline_spike_counts = bootstrapSpikeCounts(boot_data);
    
    mean_baseline_count = mean(baseline_counts);
    
    for condition = 1:numel(fits)
        x_data_stim = array_data.binEdges{condition}(stim_window_idx(1):stim_window_idx(2))'/1000; % x_data in seconds
        y_data_stim = array_data.binCounts{condition}(stim_window_idx(1)+1:stim_window_idx(2)+1)';
        
        % get response to each pulse and store
        response_to_each_pulse{condition} = zeros(numel(array_data.PULSE_TIMES{condition}{1}),1); % FR
        
        for p = 1:numel(array_data.PULSE_TIMES{condition}{1})
            spike_window_adj = spike_window/1000 + array_data.PULSE_TIMES{condition}{1}(p);
            spike_mask = array_data.spikeTrialTimes{condition} >= spike_window_adj(1) & array_data.spikeTrialTimes{condition} <= spike_window_adj(2);
            num_spikes = sum(spike_mask);
            response_to_each_pulse{condition}(p) = num_spikes/array_data.numStims(condition); % prob spike
            lat_of_spikes{condition} = [lat_of_spikes{condition}, array_data.spikeTrialTimes{condition}(spike_mask)-array_data.PULSE_TIMES{condition}{1}(p)];
        end
            
        % check if responsive
        n_pulses = 20;
        p_base = mean_baseline_count*diff(spike_window)/mode(diff(array_data.binEdges{condition})); n_base = array_data.numStims(condition)*numel(fits);
        p_stim = mean(response_to_each_pulse{condition}(1:n_pulses)); n_stim = n_pulses*array_data.numStims(condition);
        p_base = min(max(p_base,0),1);
        p_stim = min(max(p_stim,0),1);
        is_responsive(condition) = pearsonChiSquareBinomial(p_base,p_stim,n_base,n_stim,'lower') < 0.05/numel(fits); 
%             is_responsive(condition) = mean(y_data_stim(1:4)) > mean(y_data_baseline)+2*std(y_data_baseline) & mean(y_data_baseline) > min_baseline_counts;

        is_responsive_nonstim(condition ) = sum(baseline_spike_counts > p_stim)/numel(baseline_spike_counts) < 0.05/numel(fits);
        
        
        % compute response amp
        if(response_amp_num_pulses > 0) % count spikes after a certain number of pulses
            response_amp(condition) = mean(response_to_each_pulse{condition}(1:response_amp_num_pulses));
        else % count spikes in a time window, ignore spikes too close to stim onset
            num_spikes = 0;
            num_stim_pulses = 0;
            for p = 1:numel(array_data.PULSE_TIMES{condition}{1})-1
                % this one counts all spikes within in the time window
                % (except for those occuring near stim onset)
%                 spike_window_adj = [spike_window(1)/1000 + array_data.PULSE_TIMES{condition}{1}(p),array_data.PULSE_TIMES{condition}{1}(p+1)];
                % this one counts spikes within a time window post stim for 
                % each pulse in the response amp time window
                spike_window_adj = [spike_window(1)/1000 + array_data.PULSE_TIMES{condition}{1}(p),...
                    spike_window(2)/1000 + array_data.PULSE_TIMES{condition}{1}(p)];
                spike_window_adj = [min(spike_window_adj(1),response_amp_time/1000),min(spike_window_adj(2),response_amp_time/1000)];
                
                num_spikes = num_spikes + sum(array_data.spikeTrialTimes{condition} >= spike_window_adj(1) & array_data.spikeTrialTimes{condition} <= spike_window_adj(2));
                
                if(array_data.PULSE_TIMES{condition}{1}(p) < response_amp_time/1000)
                    num_stim_pulses = num_stim_pulses + 1;
                end
            end
%             response_amp(condition) = (num_spikes/array_data.numStims(condition))/(response_amp_time - num_stim_pulses*spike_window(1))*1000 - ...
%                 mean(y_data_baseline)/mode(diff(array_data.binEdges{condition})/1000); % Hz above baseline
            response_amp(condition) = (num_spikes/array_data.numStims(condition))/(num_stim_pulses*diff(spike_window))*1000 - ...
                mean(y_data_baseline)/mode(diff(array_data.binEdges{condition})/1000); % Hz above baseline

        end
        
        
        % if responsive, fit data
        if(is_responsive(condition))
            if(is_intermittent_stim==1)
                % remove bins with one or less stims in them
                keep_mask = histcounts(array_data.PULSE_TIMES{condition}{1},array_data.binEdges{condition}(stim_window_idx(1):stim_window_idx(2))/1000) > 1;
                [fits{condition},gof{condition}] = fit(x_data_stim(keep_mask),y_data_stim(keep_mask),'a*exp(-b*x)','StartPoint',[10,1],'Lower',[-1,-1]);
%                 [fits{condition},gof{condition}] = fit(x_data_stim(keep_mask),y_data_stim(keep_mask),'1+b*x','StartPoint',[0]);
            else
                [fits{condition},gof{condition}] = fit(x_data_stim,y_data_stim,'a*exp(-b*x)','StartPoint',[10,1],'Lower',[-1,-1]);
%                 [fits{condition},gof{condition}] = fit(x_data_stim(keep_mask),y_data_stim(keep_mask),'1+b*x','StartPoint',[0]);
            end
            
            param_list(condition,:) = [fits{condition}.a, fits{condition}.b];
        else
            fits{condition} = nan; gof{condition} = nan; param_list(condition,:) = [nan,nan];
        end
        
    end
    
    stim_chan = array_data.CHAN_LIST{1}(1);
    stim_chan_idx = find(map_data.chan == stim_chan);
    rec_chan_idx = find(map_data.chan == array_data.CHAN_REC);
    
    stim_row = 11 - map_data.row(stim_chan_idx);
    stim_col = map_data.col(stim_chan_idx);
    rec_row = 11 - map_data.row(rec_chan_idx);
    rec_col = map_data.col(rec_chan_idx);
    
    distance_from_stim_chan = 400*sqrt((rec_row-stim_row)^2 + (rec_col-stim_col)^2);
    
    % format outputs
    output_data.fits = fits;
    output_data.gof = gof;
    output_data.param_list = param_list;
    output_data.is_responsive = is_responsive;
    output_data.is_responsive_nonstim = is_responsive_nonstim;
    output_data.baseline_counts = baseline_counts;
    output_data.distance_from_stim_chan = distance_from_stim_chan;
    output_data.response_amp = response_amp;
    output_data.response_to_each_pulse = response_to_each_pulse;
    output_data.spike_lat = lat_of_spikes;
end

