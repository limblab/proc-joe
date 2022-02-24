function [output_data] = getReboundExcitation(array_data,input_data)

    % identify rebound excitation based on firing rate after the last pulse
    % being greater than the baseline FR by some amount

    % output whether each cell for each condition has rebound excitation
    % (1), does not (0), or did not have the condition (nan). Output
    % duration and magnitude relative to baseline FR when approapriate.
    
    % rebin (or bin for the first time)
    array_data_rebin = array_data;
    bin_edges = input_data.bin_window(1):input_data.bin_size:input_data.bin_window(2); % in ms
    for i_amp = 1:numel(array_data_rebin.binCounts)
        spike_trial_times = array_data_rebin.spikeTrialTimes{i_amp}; % in s 
        array_data_rebin.binEdges{i_amp} = bin_edges;
        array_data_rebin.binCounts{i_amp} = histcounts(spike_trial_times,bin_edges);
    end
    
    is_rebound = 0;
    rebound_dur = nan;
    rebound_mag = nan;
    
    % get threshold based on pre_window over all conditions (sets same
    % threshold for each condition)
    threshold_all = 0;
    for i_cond = input_data.cond_list
        PSTH = array_data_rebin.binCounts{i_cond};
        % smooth with running average of N kernel_length bins
        filtered_PSTH = movmean(PSTH,[input_data.kernel_length,0]);

        % compute spontaneous firing rate
        pre_window_idx = [find(array_data_rebin.binEdges{i_cond} >= input_data.pre_window(1),1,'first'),...
            find(array_data_rebin.binEdges{i_cond} >= input_data.pre_window(2),1,'first')];
        spont_fr = mean(filtered_PSTH(pre_window_idx(1):pre_window_idx(2)));

        % set threshold based on spont_fr
        threshold_all = threshold_all + spont_fr*2;
    end
    cond_mask = array_data_rebin.numStims(input_data.cond_list) > 0;
    threshold = threshold_all/sum(cond_mask); % take mean across conditions
    % get inhibition duration for this condition
    
    PSTH = array_data_rebin.binCounts{input_data.cond_idx};

    % blank during stimulation
    blank_idx = [find(array_data_rebin.binEdges{input_data.cond_idx} >= input_data.blank_time(1),1,'first'),...
        find(array_data_rebin.binEdges{input_data.cond_idx} >= input_data.blank_time(2),1,'first')];
    
    PSTH(blank_idx(1):blank_idx(2)) = mean(PSTH(pre_window_idx(1):pre_window_idx(2)));
    
    % smooth with running average of N kernel_length bins
    filtered_PSTH = movmean(PSTH,[input_data.kernel_length,0]);

    % if firing rate overshoots thresh, define duration as the time
    % between this point and when the FR comes back above thresh. 
    post_window_idx = [find(array_data_rebin.binEdges{input_data.cond_idx} >= input_data.post_window(1),1,'first'),...
        find(array_data_rebin.binEdges{input_data.cond_idx} >= input_data.post_window(2),1,'first')];
    over_thresh_mask(1,:) = squeeze(filtered_PSTH(post_window_idx(1):post_window_idx(2)) > threshold);

    start_ones = strfind([over_thresh_mask,1],[0 1]);
    end_ones = strfind([over_thresh_mask==1,0],[1 0]);
    if(~isempty(end_ones) && ~isempty(start_ones) && end_ones(1) < start_ones(1))
        end_ones(1) = [];
    end
    num_consecutive_ones = end_ones(1:min(numel(end_ones),numel(start_ones))) - start_ones(1:min(numel(end_ones),numel(start_ones)));

    [num_consecutive_ones,max_idx] = max(num_consecutive_ones);
    if(~isempty(num_consecutive_ones) && num_consecutive_ones >= input_data.num_consec_bins)
        is_rebound = 1;
        rebound_dur = num_consecutive_ones*input_data.bin_size;
        rebound_mag = 0;
    end      
    
    output_data = [];
    output_data.PSTH = PSTH;
    output_data.filtered_PSTH = filtered_PSTH;
    output_data.is_rebound = is_rebound;
    output_data.rebound_dur = rebound_dur;
    output_data.rebound_mag = rebound_mag;
    output_data.threshold = threshold;
end