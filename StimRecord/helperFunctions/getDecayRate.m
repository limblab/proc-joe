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
    baseline_counts = zeros(numel(fits),1);
    response_to_each_pulse = cell(numel(fits),1);
    
    response_amp = zeros(numel(fits),1);
    
    for condition = 1:numel(fits)
        x_data_stim = array_data.binEdges{condition}(stim_window_idx(1):stim_window_idx(2))'/1000; % x_data in seconds
        y_data_stim = array_data.binCounts{condition}(stim_window_idx(1):stim_window_idx(2))';
        
        x_data_baseline = array_data.binEdges{condition}(baseline_window_idx(1):baseline_window_idx(2))'/1000; % in seconds
        y_data_baseline = array_data.binCounts{condition}(baseline_window_idx(1):baseline_window_idx(2))';
        
        % check if responsive
        is_responsive(condition) = mean(y_data_stim(1:4)) > mean(y_data_baseline)+2*std(y_data_baseline) & mean(y_data_baseline) > min_baseline_counts;
        baseline_counts(condition) = mean(y_data_baseline);
        
        % get response to each pulse and store
        response_to_each_pulse{condition} = zeros(numel(array_data.PULSE_TIMES{condition}{1}),1); % FR
        
        for p = 1:numel(array_data.PULSE_TIMES{condition}{1})
            spike_window_adj = spike_window/1000 + array_data.PULSE_TIMES{condition}{1}(p);
            num_spikes = sum(array_data.spikeTrialTimes{condition} >= spike_window_adj(1) & array_data.spikeTrialTimes{condition} <= spike_window_adj(2));
            response_to_each_pulse{condition}(p) = num_spikes/(diff(spike_window/1000))/array_data.numStims(condition) - ...
                mean(y_data_baseline)/mode(diff(array_data.binEdges{condition})/1000); % Hz above baseline
        end
            
        % compute response amp
        if(response_amp_num_pulses > 0) % count spikes after a certain number of pulses
            response_amp(condition) = mean(response_to_each_pulse{condition}(1:response_amp_num_pulses));
        else % count spikes in a time window, ignore spikes too close to stim onset
            num_spikes = 0;
            num_stim_pulses = 0;
            for p = 1:numel(array_data.PULSE_TIMES{condition}{1})-1
                spike_window_adj = [spike_window(1)/1000 + array_data.PULSE_TIMES{condition}{1}(p),array_data.PULSE_TIMES{condition}{1}(p+1)];
                spike_window_adj = [min(spike_window_adj(1),response_amp_time/1000),min(spike_window_adj(2),response_amp_time/1000)];
                
                num_spikes = num_spikes + sum(array_data.spikeTrialTimes{condition} >= spike_window_adj(1) & array_data.spikeTrialTimes{condition} <= spike_window_adj(2));
                
                if(array_data.PULSE_TIMES{condition}{1}(p) < response_amp_time/1000)
                    num_stim_pulses = num_stim_pulses + 1;
                end
            end
            response_amp(condition) = (num_spikes/array_data.numStims(condition))/(response_amp_time - num_stim_pulses*spike_window(1))*1000 - ...
                mean(y_data_baseline)/mode(diff(array_data.binEdges{condition})/1000); % Hz above baseline
        end
        
        
        % if responsive, fit data
        if(is_responsive(condition))
            if(is_intermittent_stim==1)
                % remove bins with one or less stims in them
                keep_mask = histcounts(array_data.PULSE_TIMES{condition}{1},array_data.binEdges{condition}(stim_window_idx(1):stim_window_idx(2))/1000) > 1;
                [fits{condition},gof{condition}] = fit(x_data_stim(keep_mask),y_data_stim(keep_mask),'a*exp(-b*x)','StartPoint',[10,1],'Lower',[-10,-10]);
            else
                [fits{condition},gof{condition}] = fit(x_data_stim,y_data_stim,'a*exp(-b*x)','StartPoint',[10,1],'Lower',[0,0]);
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
    output_data.baseline_counts = baseline_counts;
    output_data.distance_from_stim_chan = distance_from_stim_chan;
    output_data.response_amp = response_amp;
    output_data.response_to_each_pulse = response_to_each_pulse;
end

