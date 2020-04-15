function [output_data] = getReboundExcitationStats(array_data,input_data)

% RETURNS rebound_amp, peak_time, onset_latency, duration if a rebound occurs, otherwise
% these are empty -- does this for each condition separately
% rebound_amp = integrate above baseline for the duration
% peak_time = peak if it exists
% duration = time from onset to offset
    
    rebound_amp = nan(size(array_data.spikeTrialTimes));
    peak_time = nan(size(rebound_amp));
    duration = nan(size(rebound_amp));
    is_rebound = zeros(size(rebound_amp));
    
    % rebin array data according to bin size
    
    array_data = rebinArrayData(array_data,input_data.bin_size);
    
    % get mean and std baseline fr across all conditions
    bin_edge_idx = [find(array_data.binEdges{1} > input_data.baseline_window(1),1,'first'),...
        find(array_data.binEdges{1} > input_data.baseline_window(2),1,'first')];
    
    baseline_counts = [];
    
    for condition = 1:numel(array_data.spikeTrialTimes)
        baseline_counts = [baseline_counts, array_data.binCounts{condition}(bin_edge_idx(1):bin_edge_idx(2))];
    end
    mean_baseline_count = mean(baseline_counts,'omitnan');
    std_baseline_count = std(baseline_counts,'omitnan');
    
    % for each condition
    for condition = 1:numel(array_data.spikeTrialTimes)
        if(array_data.num_stims(condition) > 0) % do stuff
            % determine if there is a rebound excitation
            % find peak in window provided, if peak is bigger than mean +
            % thresh*std, we have rebound excitation

            % fix a bug where the pulse time is longer than the actual
            % train for some reason
            pulse_time_idx = 1;
            while(pulse_time_idx <= numel(array_data.PULSE_TIMES{condition}) && array_data.PULSE_TIMES{condition}{pulse_time_idx}(end)*1000 > 250)
                pulse_time_idx = pulse_time_idx + 1;
            end
            
            if(pulse_time_idx <= numel(array_data.PULSE_TIMES{condition}))
                post_stim_idx = [find(array_data.binEdges{condition} > array_data.PULSE_TIMES{condition}{pulse_time_idx}(end)*1000 + input_data.post_stim_window(1),1,'first'),...
                    find(array_data.binEdges{condition} > array_data.PULSE_TIMES{condition}{pulse_time_idx}(end)*1000 + input_data.post_stim_window(2),1,'first')];

                above_threshold = array_data.binCounts{condition}(post_stim_idx(1):post_stim_idx(2)) > mean_baseline_count + input_data.threshold_mult*std_baseline_count;

                up_step_idx = strfind([0,above_threshold],[0 1]);
                consecutive_ones = strfind([above_threshold,0],[1 0]) - up_step_idx + 1;

                [max_consec, max_idx] = max(consecutive_ones);

                if(~isempty(max_consec))
                    is_rebound(condition) = max_consec  > input_data.num_bins_above_thresh;
                end

                if(is_rebound(condition)) % compute relevant statistics
                    [~, peak_idx] = max(array_data.binCounts{condition}(post_stim_idx(1):post_stim_idx(2)));
                    zero_idx = find(array_data.binEdges{condition} > 0,1,'first');

                    peak_time(condition) = (peak_idx - 1)*input_data.bin_size + input_data.post_stim_window(1);
                    duration(condition) = max_consec*input_data.bin_size;
                    rebound_amp(condition) = sum(array_data.binCounts{condition});

                end
            end
        end
    end
    
    
    output_data.rebound_amp = rebound_amp;
    output_data.peak_time = peak_time;
    output_data.duration = duration;
    output_data.is_rebound = is_rebound;

end