function [array_data] = getResponseAmplitude(array_data, bin_size)

    array_data = rebinArrayData(array_data,bin_size);

    is_responsive = zeros(size(array_data.binCounts));
    response_amp = zeros(size(array_data.binCounts));
    
    for chan = 1:size(array_data.binCounts,1)
        for wave = 1:size(array_data.binCounts,2)
            % find peak period (if it exists) -- 3 consecutive bins above mean baseline +
            % 2*SD baseline

            baseline_bins = [find(array_data.binEdges{chan,wave} > -20,1,'first'), find(array_data.binEdges{chan,wave} > -2,1,'first')];
            zero_bin = find(array_data.binEdges{chan,wave} > 0,1,'first');
            stim_window_idx = [find(array_data.binEdges{chan,wave} > 1,1,'first'),...
                find(array_data.binEdges{chan,wave} > 10,1,'first')];

            mean_baseline = mean(array_data.binCounts{chan,wave}(baseline_bins(1):baseline_bins(2)));
            std_baseline = std(array_data.binCounts{chan,wave}(baseline_bins(1):baseline_bins(2)));

            threshold = mean_baseline + 2*std_baseline;

            above_threshold = array_data.binCounts{chan,wave}(stim_window_idx(1):stim_window_idx(2)) > threshold;
            up_step_idx = strfind([0,above_threshold],[0 1]);
            consecutive_ones = strfind([above_threshold,0],[1 0]) - up_step_idx + 1;

            [max_consec, max_idx] = max(consecutive_ones);
            % if peak period, sum spikes in that region
            if(max_consec >= 3)
                response_amp(chan,wave) = (sum(array_data.binCounts{chan,wave}(stim_window_idx(1)+up_step_idx(max_idx):...
                    stim_window_idx(1)+up_step_idx(max_idx)+max_consec-1)) - mean_baseline*max_consec)/array_data.numStims(chan,wave);
                is_responsive(chan,wave) = 1;


            else % else sum spikes in 3 bin window centered at the peak
                [~,max_idx] = max(array_data.binCounts{chan,wave}(stim_window_idx(1):stim_window_idx(2)));
                
                response_amp(chan,wave) = (sum(array_data.binCounts{chan,wave}(stim_window_idx(1)+max_idx-2:stim_window_idx(1)+max_idx)) - ...
                    mean_baseline*3)/array_data.numStims(chan,wave);
                is_responsive(chan,wave) = 0;

            end
        
        end
        
    end
    
    
    array_data.is_responsive = is_responsive;
    array_data.response_amp = response_amp;

 


end