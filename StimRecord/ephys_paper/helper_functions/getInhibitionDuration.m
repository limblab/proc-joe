function [output_data] = getInhibitionDuration(array_data,cond_idx,input_data)

    % gets inhibition duration for the i_ampitions in array_data. Inhibition
    % duration is defined based on (Butovas 2003)
    
    % input_data contains
    %   PRE_WINDOW window used to compute spontaneous firing rate in ms
    %   POST_WINDOW window used to compute inhib duration (if exists)
    %   BLANK_TIME -- time after stimulation onset to blank to remove
    %       excitatory period due to stimulation
      
    filtered_PSTH = [];
    PSTH = [];
    threshold = [];
    is_inhib = 0;
    inhib_dur = nan;
    inhib_off_time = nan;
    
    % rebin (or bin for the first time)
    array_data_rebin = array_data;
    bin_edges = input_data.bin_window(1):input_data.bin_size:input_data.bin_window(2); % in ms
    for i_amp = 1:numel(array_data_rebin.binCounts)
        spike_trial_times = array_data_rebin.spikeTrialTimes{i_amp}; % in s 
        array_data_rebin.binEdges{i_amp} = bin_edges;
        array_data_rebin.binCounts{i_amp} = histcounts(spike_trial_times,bin_edges);
    end
    
    % get threshold based on pre_window over all conditions (sets same
    % threshold for each condition)
    threshold_all = 0;
    for i_cond = 1:numel(array_data_rebin.binCounts)
        PSTH = array_data_rebin.binCounts{1,i_cond};
        % smooth with running average of N kernel_length bins
        filtered_PSTH = movmean(PSTH,[input_data.kernel_length,0]);

        % compute spontaneous firing rate
        pre_window_idx = [find(array_data_rebin.binEdges{1,i_cond} >= input_data.pre_window(1),1,'first'),...
            find(array_data_rebin.binEdges{1,i_cond} >= input_data.pre_window(2),1,'first')];
        spont_fr = mean(filtered_PSTH(pre_window_idx(1):pre_window_idx(2)));
        spont_std = std(filtered_PSTH(pre_window_idx(1):pre_window_idx(2)));

        % set threshold based on spont_fr
        threshold_all = threshold_all + spont_fr*0.75;
    end
    cond_mask = array_data_rebin.numStims > 0;
    threshold = threshold_all/sum(cond_mask); % take mean across conditions
    % get inhibition duration for this condition
    
    PSTH = array_data_rebin.binCounts{1,cond_idx};

    % blank period with mean from baseline
    blank_idx = [find(array_data_rebin.binEdges{1,cond_idx} >= input_data.blank_time(1),1,'first'),...
        find(array_data_rebin.binEdges{1,cond_idx} >= input_data.blank_time(2),1,'first')];
    pre_window_idx = [find(array_data_rebin.binEdges{1,cond_idx} >= input_data.pre_window(1),1,'first'),...
        find(array_data_rebin.binEdges{1,cond_idx} >= input_data.pre_window(2),1,'first')];

    PSTH(blank_idx(1):blank_idx(2)) = mean(PSTH(pre_window_idx(1):pre_window_idx(2)));

%         % filter with gaussian kernel (for Butovas case)
%         filtered_PSTH(i_amp,:) = gaussianKernel(PSTH(i_amp,:),input_data.kernel_length);

    % smooth with running average of N kernel_length bins
    filtered_PSTH = movmean(PSTH,[input_data.kernel_length,0]);

     % if firing rate undershoots thresh, define duration as the time
    % between this point and when the FR comes back above thresh. 
    post_window_idx = [find(array_data_rebin.binEdges{1,cond_idx} >= input_data.post_window(1),1,'first'),...
        find(array_data_rebin.binEdges{1,cond_idx} >= input_data.post_window(2),1,'first')];
    under_thresh_mask(1,:) = squeeze(filtered_PSTH(post_window_idx(1):post_window_idx(2)) < threshold);

    start_ones = strfind([under_thresh_mask,1],[0 1]);
    start_ones(start_ones > input_data.max_time_start/input_data.bin_size) = [];
    end_ones = strfind([under_thresh_mask==1,0],[1 0]);
    if(~isempty(end_ones) && ~isempty(start_ones) && end_ones(1) < start_ones(1))
        end_ones(1) = [];
    end
    num_consecutive_ones = end_ones(1:min(numel(end_ones),numel(start_ones))) - start_ones(1:min(numel(end_ones),numel(start_ones))) + 1;

    [num_consecutive_ones,max_idx] = max(num_consecutive_ones);
    if(~isempty(num_consecutive_ones) && num_consecutive_ones >= input_data.num_consec_bins)
        is_inhib = 1;
        inhib_dur = num_consecutive_ones*input_data.bin_size;
        inhib_off_time = end_ones(max_idx)*input_data.bin_size;
    end      
    
    output_data.filtered_PSTH = filtered_PSTH;
    output_data.is_inhib = is_inhib;
    output_data.inhib_dur = inhib_dur;
    output_data.PSTH = PSTH;
    output_data.threshold = threshold;
    output_data.inhib_off_time = inhib_off_time;
end



function data_smooth = gaussianKernel(data,kernel_length)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nbr_samples = numel(data);
    % apply smoothing to data

    % kernel half length is 3·SD out
    kernel_hl = ceil( 3 * kernel_length  );
    % create the kernel --it will have length 2*kernel_hl+1
    kernel = normpdf( -kernel_hl : ...
        1 : kernel_hl, ...
        0, kernel_length );
    
    % compute normalization factor --this factor depends on the number of taps
    % actually used
    nm = conv(kernel,ones(1,nbr_samples));

    % do the smoothing
        aux_smoothed_FR     = conv(kernel,data) ./ nm;
        % cut off the edges so that the result of conv is same length as the
        % original data
        data_smooth    = aux_smoothed_FR(kernel_hl+1:end-kernel_hl);
end
