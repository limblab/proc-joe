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
    amp = [];
    threshold = [];
    is_inhib = zeros(numel(array_data_rebin.binCounts),1);
    inhib_dur = nan(numel(array_data_rebin.binCounts),1);
    % for each condition:
    for cond = 1:numel(array_data_rebin.binCounts)
    % low pass filter PSTH with a gaussian kernal, length 50 bins (1 ms bin
    % width)
        PSTH(cond,:) = array_data_rebin.binCounts{cond};
        
        % blank period with mean from baseline
        blank_idx = [find(array_data_rebin.binEdges{cond} >= input_data.BLANK_TIME(1),1,'first'),...
            find(array_data_rebin.binEdges{cond} > input_data.BLANK_TIME(2),1,'first')];
        pre_window_idx = [find(array_data_rebin.binEdges{cond} > input_data.PRE_WINDOW(1),1,'first'),...
            find(array_data_rebin.binEdges{cond} > input_data.PRE_WINDOW(2),1,'first')];
        
        PSTH(cond,blank_idx(1):blank_idx(2)) = mean(PSTH(cond,pre_window_idx(1):pre_window_idx(2)));
        
%         % filter with gaussian kernel (for Butovas case)
%         filtered_PSTH(cond,:) = gaussianKernel(array_data_rebin.binCounts{cond},input_data.KERNEL_LENGTH);
        
        % smooth with running average of N kernel_length bins
        filtered_PSTH(cond,:) = movmean(PSTH(cond,:),[input_data.KERNEL_LENGTH,0]);
        
%         % compute spontaneous firing rate
%         spont_fr = mean(filtered_PSTH(cond,pre_window_idx(1):pre_window_idx(2)));
%         spont_std = std(filtered_PSTH(cond,pre_window_idx(1):pre_window_idx(2)));
%         
%         % set threshold based on standard deviation
%         threshold(cond) = spont_fr - spont_std;
% %         threshold(cond) = spont_fr*0.75;

    end
    
    for cond = 1:numel(array_data_rebin.binCounts)
        spont_fr = mean(reshape(filtered_PSTH(cond,pre_window_idx(1):pre_window_idx(2)),numel(filtered_PSTH(cond,pre_window_idx(1):pre_window_idx(2))),1));
        spont_std = std(reshape(filtered_PSTH(cond,pre_window_idx(1):pre_window_idx(2)),numel(filtered_PSTH(cond,pre_window_idx(1):pre_window_idx(2))),1));
        
        
        threshold = zeros(numel(array_data_rebin.binCounts),1)+spont_fr - spont_std;
        
        % if firing rate undershoots thresh, define duration as the time
        % between this point and when the FR comes back above thresh. 
        post_window_idx = [find(array_data_rebin.binEdges{cond} > input_data.POST_WINDOW(1),1,'first'),...
            find(array_data_rebin.binEdges{cond} > input_data.POST_WINDOW(2),1,'first')];
        under_one_mask = (filtered_PSTH(cond,post_window_idx(1):post_window_idx(2)) < threshold(cond));
        
        start_ones = strfind([under_one_mask,1],[0 1]);
        start_ones(start_ones > input_data.MAX_TIME_START) = [];
        end_ones = strfind([under_one_mask==1,0],[1 0]);
        if(~isempty(end_ones) && ~isempty(start_ones) && end_ones(1) < start_ones(1))
            end_ones(1) = [];
        end
        num_consecutive_ones = end_ones(1:min(numel(end_ones),numel(start_ones))) - start_ones(1:min(numel(end_ones),numel(start_ones))) + 1;
        
        num_consecutive_ones = max(num_consecutive_ones);
        if(~isempty(num_consecutive_ones) && num_consecutive_ones >= input_data.NUM_CONSECUTIVE_BINS)
            is_inhib(cond) = 1;
            inhib_dur(cond) = num_consecutive_ones*input_data.BIN_SIZE;
        end
%         if(~isempty(under_thresh_idx) && under_thresh_idx(1)+post_window_idx(1) < find(array_data_rebin.binEdges{cond} > input_data.MAX_TIME_START,1,'first'))
%             is_inhib(cond) = 1;
%             
%             % find end of first chain
%             chain_idx = find(diff(under_thresh_idx) > 1,1,'first');
%             if(isempty(chain_idx))
%                 inhib_dur(cond) = (under_thresh_idx(end)-under_thresh_idx(1))*input_data.BIN_SIZE;
%             else
%                 inhib_dur(cond) = (under_thresh_idx(chain_idx)-under_thresh_idx(1))*input_data.BIN_SIZE;
%             end
%         end
        
        
    end
    
    % check against parameters
    if(isfield(array_data,'STIM_PARAMETERS'))
        keep_mask = zeros(size(inhib_dur)); % determines which conditions to keep
        for cond = 1:numel(array_data.binCounts)
            amp(cond) = array_data.STIM_PARAMETERS(cond).amp1;

            if(array_data.STIM_PARAMETERS(cond).polarity == input_data.POL && ...
                    array_data.STIM_PARAMETERS(cond).pWidth1 == input_data.PW1 && ...
                    array_data.STIM_PARAMETERS(cond).pWidth2 == input_data.PW2 && ...
                    is_inhib(cond)==1)
                % then keep this condition
                keep_mask(cond) = 1;
            end

        end
    else
        keep_mask = []; amp = [];
    end
%     plot(amp(keep_mask==1),inhib_dur(keep_mask==1),'marker','.','markersize',16)
    
    inhib_struct.filtered_PSTH = filtered_PSTH;
    inhib_struct.is_inhib = is_inhib;
    inhib_struct.inhib_dur = inhib_dur;
    inhib_struct.PSTH = PSTH;
    inhib_struct.threshold = threshold;
    inhib_struct.amp = amp;
    inhib_struct.keep_mask = keep_mask;
    
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
