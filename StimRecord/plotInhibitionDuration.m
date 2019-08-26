function [inhib_struct, figure_handles] = plotInhibitionDuration(array_data,input_data)

    % gets inhibition duration for the conditions in array_data. Inhibition
    % duration is defined based on (Butovas 2003)
    
    % input_data contains
    %   PRE_WINDOW window used to compute spontaneous firing rate in ms
    %   POST_WINDOW window used to compute inhib duration (if exists)
    %   BLANK_TIME -- time after stimulation onset to blank to remove
    %       excitatory period due to stimulation
    
    
    figure_handles = [];
    
    % rebin to 1 ms bins
    array_data_rebin = array_data;
    bin_edges = array_data_rebin.binEdges{1}(1):input_data.BIN_SIZE:array_data_rebin.binEdges{1}(end); % in ms
    for cond = 1:numel(array_data_rebin.binCounts)
        spike_trial_times = array_data_rebin.spikeTrialTimes{cond}; % in s for some reason
        array_data_rebin.binEdges{cond} = bin_edges;
        array_data_rebin.binCounts{cond} = histcounts(spike_trial_times*1000,bin_edges);
    end

    
    filtered_PSTH = [];
    PSTH = [];
    spont_fr = [];
    is_inhib = zeros(numel(array_data_rebin.binCounts),1);
    inhib_dur = -1+zeros(numel(array_data_rebin.binCounts),1);
    % for each condition:
    for cond = 1:numel(array_data_rebin.binCounts)
    % low pass filter PSTH with a gaussian kernal, length 50 bins (1 ms bin
    % width)
        PSTH(cond,:) = array_data_rebin.binCounts{cond};
        
        % blank period with mean from baseline
        blank_idx = [find(array_data_rebin.binEdges{cond} >= 0,1,'first'),...
            find(array_data_rebin.binEdges{cond} > input_data.BLANK_TIME,1,'first')];
        pre_window_idx = [find(array_data_rebin.binEdges{cond} > input_data.PRE_WINDOW(1),1,'first'),...
            find(array_data_rebin.binEdges{cond} > input_data.PRE_WINDOW(2),1,'first')];
        
        PSTH(cond,blank_idx(1):blank_idx(2)) = mean(PSTH(cond,pre_window_idx(1):pre_window_idx(2)));
        
        % filter with gaussian kernel
        filtered_PSTH(cond,:) = gaussianKernel(array_data_rebin.binCounts{cond},input_data.KERNEL_LENGTH);
        
        % compute spontaneous firing rate
        spont_fr(cond) = mean(filtered_PSTH(cond,pre_window_idx(1):pre_window_idx(2)));
        
    end
    % set threshold as 0.75 mean(fr_spontaneous)
    threshold = 0.75*mean(spont_fr);
    
    for cond = 1:numel(array_data_rebin.binCounts)
        % if firing rate undershoots thresh, define duration as the time
        % between this point and when the FR comes back above thresh. 
        post_window_idx = [find(array_data_rebin.binEdges{cond} > input_data.POST_WINDOW(1),1,'first'),...
            find(array_data_rebin.binEdges{cond} > input_data.POST_WINDOW(2),1,'first')];
        under_thresh_idx = find(filtered_PSTH(cond,post_window_idx(1):post_window_idx(2)) < threshold);
        
        if(~isempty(under_thresh_idx))
            is_inhib(cond) = 1;
            
            % find end of first chain
            chain_idx = find(diff(under_thresh_idx) > 1);
            if(isempty(chain_idx))
                inhib_dur(cond) = under_thresh_idx(end)*input_data.BIN_SIZE;
            else
                inhib_dur(cond) = under_thresh_idx(chain_idx)*input_data.BIN_SIZE;
            end
        end
        
    end
    
    
    
    inhib_struct.filtered_PSTH = filtered_PSTH;
    inhib_struct.is_inhib = is_inhib;
    inhib_struct.inhib_dur = inhib_dur;
    inhib_struct.PSTH = PSTH;
    inhib_struct.threshold = threshold;

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
