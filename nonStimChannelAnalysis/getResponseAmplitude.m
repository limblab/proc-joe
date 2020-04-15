function [array_data] = getResponseAmplitude(array_data, arr_idx,bin_size,window_width)

    array_data = rebinArrayData(array_data,bin_size);

    is_responsive = zeros(size(array_data.binCounts));
    response_amp = zeros(size(array_data.binCounts));
    peak_latency = zeros(size(array_data.binCounts));
    mean_baseline = zeros(size(array_data.binCounts));
    std_baseline = zeros(size(array_data.binCounts));
    
    for chan = 1:size(array_data.binCounts,1)
        for wave = 1:size(array_data.binCounts,2)
            % find peak period (if it exists) -- 3 consecutive bins above mean baseline +
            % 2*SD baseline

            baseline_bins = [find(array_data.binEdges{chan,wave} > -20,1,'first'), find(array_data.binEdges{chan,wave} > -2,1,'first')];
            zero_bin = find(array_data.binEdges{chan,wave} > 0,1,'first');
            stim_window_idx = [find(array_data.binEdges{chan,wave} > 1,1,'first'),...
                find(array_data.binEdges{chan,wave} > 7,1,'first')];

            mean_baseline(chan,wave) = mean(array_data.binCounts{chan,wave}(baseline_bins(1):baseline_bins(2)));
            std_baseline(chan,wave) = std(array_data.binCounts{chan,wave}(baseline_bins(1):baseline_bins(2)));

            threshold = mean_baseline(chan,wave) + 2*std_baseline(chan,wave);

            above_threshold = array_data.binCounts{chan,wave}(stim_window_idx(1):stim_window_idx(2)) > threshold;
            up_step_idx = strfind([0,above_threshold],[0 1]);
            consecutive_ones = strfind([above_threshold,0],[1 0]) - up_step_idx + 1;

            [max_consec, max_idx] = max(consecutive_ones);
            % if peak period, sum spikes in that region
            if(~isempty(max_consec) && max_consec >= 3 && mean_baseline(chan,wave) > 0.001)
                response_amp(chan,wave) = (sum(array_data.binCounts{chan,wave}(stim_window_idx(1)+up_step_idx(max_idx)-1:...
                    stim_window_idx(1)+up_step_idx(max_idx)+max_consec-2)) - mean_baseline(chan,wave)*max_consec);
                is_responsive(chan,wave) = 1;
                
                % get peak latency 
                max_idx = find(consecutive_ones >= 3,1,'first');
                
                [~,first_peak_idx] = max(array_data.binCounts{chan,wave}(stim_window_idx(1)+up_step_idx(max_idx)-1:stim_window_idx(1)+up_step_idx(max_idx)+consecutive_ones(1)));
%                 [~,first_peak_idx] = findpeaks(array_data.binCounts{chan,wave}(stim_window_idx(1)+up_step_idx(max_idx)-2:stim_window_idx(1)+up_step_idx(max_idx)+consecutive_ones(1)),'NPeaks',1);
                
                peak_latency(chan,wave) = (up_step_idx(max_idx) - 1 + first_peak_idx - 1)*bin_size + 1 - 0.453; % end of stimulation
                
                if(peak_latency(chan,wave) > 3 && wave == 4 && array_data.mean_time_off_rail(chan,wave) <= 1)
%                    figure();
%                    plot(array_data.binCounts{chan,wave}(stim_window_idx(1):stim_window_idx(2)))
%                     disp(mean_baseline(chan,wave))
%                     disp([num2str(arr_idx),' ',num2str(chan),' ',num2str(wave),' ',num2str(peak_latency(chan,wave))])
                end
  
            else
                
%             sum spikes in 3 bin window centered at the peak
                [~,max_idx] = max(array_data.binCounts{chan,wave}(stim_window_idx(1):stim_window_idx(2)));

%                 response_amp(chan,wave) = (sum(array_data.binCounts{chan,wave}(stim_window_idx(1)+max_idx-window_width-1:stim_window_idx(1)+max_idx+window_width-1)) - ...
%                     mean_baseline(chan,wave)*(1+2*window_width));
                response_amp(chan,wave) = (sum(array_data.binCounts{chan,wave}(stim_window_idx(1):stim_window_idx(2))) - ...
                    mean_baseline(chan,wave)*(1+diff(stim_window_idx)));
                
                is_responsive(chan,wave) = 0;
                peak_latency(chan,wave) = nan;
            end
            
              
            % sanitize data
           if(peak_latency(chan,wave) > 6 || (isfield(array_data,'mean_time_off_rail') && array_data.mean_time_off_rail(chan,wave) > 1))
                peak_latency(chan,wave) = nan;
                response_amp(chan,wave) = nan;
                is_responsive(chan,wave) = 0;
                disp(arr_idx)
           end 
%             end
        
        end
        
    end
    
    
    array_data.is_responsive = is_responsive;
    array_data.response_amp = response_amp;
    array_data.peak_latency = peak_latency;
    array_data.mean_baseline = mean_baseline;


end