function [output_data,figure_handles] = plotLatencyActivation(array_data, input_data)


    figure_handles = [];
    
    % rebin to bin size
    array_data_rebin = array_data;
    bin_edges = array_data_rebin.binEdges{1}(1):input_data.BIN_SIZE:array_data_rebin.binEdges{1}(end); % in ms
    for cond = 1:numel(array_data_rebin.binCounts)
        spike_trial_times = array_data_rebin.spikeTrialTimes{cond}; % in s for some reason
        array_data_rebin.binEdges{cond} = bin_edges;
        array_data_rebin.binCounts{cond} = histcounts(spike_trial_times*1000,bin_edges)/array_data_rebin.numStims(cond);
    end

    
    filtered_PSTH = [];
    PSTH = [];
    amp = [];
    latencies = {};
    threshold = [];
    % for each condition:
    for cond = 1:numel(array_data_rebin.binCounts)
    % low pass filter PSTH with a gaussian kernal, length 50 bins (1 ms bin
    % width)
    
        PSTH(cond,:) = array_data_rebin.binCounts{cond};
                
        % filter with gaussian kernel 
        filtered_PSTH(cond,:) = gaussianKernel(array_data_rebin.binCounts{cond},input_data.KERNEL_LENGTH);
        
        % get peak latencies in filtered_PSTH
        peak_window_idx = [find(array_data_rebin.binEdges{cond} > input_data.PEAK_WINDOW(1),1,'first'),...
            find(array_data_rebin.binEdges{cond} > input_data.PEAK_WINDOW(2),1,'first')];
        
        baseline_window_idx = [find(array_data_rebin.binEdges{cond} > input_data.BASELINE_WINDOW(1),1,'first'),...
            find(array_data_rebin.binEdges{cond} > input_data.BASELINE_WINDOW(2),1,'first')];
        
        baseline_peak_mag = findpeaks(filtered_PSTH(cond,baseline_window_idx(1):baseline_window_idx(2)));
        
        threshold(cond) = mean(baseline_peak_mag) + 3*std(baseline_peak_mag);
        
        [~,peaks] = findpeaks(filtered_PSTH(cond,peak_window_idx(1):peak_window_idx(2)),'MinPeakHeight',threshold(cond),'MinPeakDistance',ceil(1/input_data.BIN_SIZE));
        latencies{cond} = peaks*input_data.BIN_SIZE;
        
        %% adjust latencies based on stim length
        if(input_data.ADJUST_LATENCY_TIME)
            latencies{cond} = latencies{cond} - (array_data_rebin.STIM_PARAMETERS(cond).pWidth1 + ...
                array_data_rebin.STIM_PARAMETERS(cond).pWidth2 + array_data_rebin.STIM_PARAMETERS(cond).interphase)*(10^(-3)); 
        end
        
    end
    
    % check against parameters
    keep_mask = zeros(size(threshold)); % determines which conditions to keep
    for cond = 1:numel(array_data.binCounts)
        amp(cond) = array_data.STIM_PARAMETERS(cond).amp1;
        pw1(cond) = array_data.STIM_PARAMETERS(cond).pWidth1;
        pol(cond) = array_data.STIM_PARAMETERS(cond).polarity;
        
        if(array_data.STIM_PARAMETERS(cond).polarity == input_data.POL && ...
                array_data.STIM_PARAMETERS(cond).pWidth1 == input_data.PW1 && ...
                array_data.STIM_PARAMETERS(cond).pWidth2 == input_data.PW2)
            % then keep this condition
            keep_mask(cond) = 1;
        end

    end
    
    output_data.filtered_PSTH = filtered_PSTH;
    output_data.PSTH = PSTH;
    output_data.amp = amp;
    output_data.latencies = latencies;
    output_data.threshold = threshold;
    output_data.keep_mask = keep_mask;
    output_data.pw1 = pw1;
    output_data.pol = pol;
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
