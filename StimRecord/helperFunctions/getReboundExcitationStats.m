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
    mean_baseline_count = mean(baseline_counts);
    std_baseline_count = std(baseline_counts);
    
    % for each condition
    for condition = 1:numel(array_data.spikeTrialTimes)
        
        % determine if there is a rebound excitation
        % find peak in window provided, if peak is bigger than mean +
        % thresh*std, we have rebound excitation
        
        post_stim_idx = [find(array_data.binEdges{condition} > input_data.train_length(condition) + input_data.post_stim_window(1),1,'first'),...
            find(array_data.binEdges{condition} > input_data.train_length(condition) + input_data.post_stim_window(2),1,'first')];
        
        above_thresh = array_data.binCounts{condition}(post_stim_idx(1):post_stim_idx(2)) > mean_baseline_count + input_data.threshold_mult*std_baseline_count;
        
        flip_on = find(diff(above_thresh) == 1);
        flip_off = find(diff(above_thresh) == -1);
        
        
        is_rebound(condition) = sum(above_thresh)  > input_data.num_bins_above_thresh;
        
        if(is_rebound(condition)) % compute relevant statistics
            [~, peak_idx] = max(array_data.binCounts{condition}(post_stim_idx(1):post_stim_idx(2)));
            zero_idx = find(array_data.binEdges{condition} > 0,1,'first');
            
            peak_time(condition) = (peak_idx - 1)*input_data.bin_size + input_data.post_stim_window(1);
            
            
            
%             rebound_amp = nan(size(array_data.spikeTrialTimes));
%             peak_time = nan(size(rebound_amp));
%             duration = nan(size(rebound_amp));
            
        end
    end
    
    
    output_data.rebound_amp = rebound_amp;
    output_data.peak_time = peak_time;
    output_data.duration = duration;
    output_data.is_rebound = is_rebound;

end