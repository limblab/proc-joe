function [ array_data_trim ] = splitArrayData( array_data,speed_bounds, window, normalizeBinCounts )
% returns an arrayData with stimulations where the mean speed in window is
% within speed_bounds

% window is in ms relative to stimulation onset. 

    if(~isfield(array_data{1}.kin{1},'speed'))
        % compute speed
        for u = 1:numel(array_data)
            for cond = 1:numel(array_data{u}.binCounts)
                
                array_data{u}.kin{cond}.speed = sqrt((array_data{u}.kin{cond}.vx).^2 + (array_data{u}.kin{cond}.vy).^2);
                % compute mean speed within window
                window_idx = [find(array_data{u}.binEdges{cond} > window(1),1,'first'), ...
                    find(array_data{u}.binEdges{cond} > window(2),1,'first')];
                array_data{u}.kin{cond}.mean_speed = mean(array_data{u}.kin{cond}.speed(:,window_idx(1):window_idx(2)),2);
            end
        end
        
    end

    array_data_trim = array_data;
    for u = 1:numel(array_data)
        for cond = 1:numel(array_data{u}.binCounts)
            
            if(isfield(array_data{u},'trial_num')) % because I called the same thing a different name for trains and single pulses....
                trial_speeds = array_data{u}.kin{cond}.mean_speed(array_data{u}.trial_num{cond});
            else
                trial_speeds = array_data{u}.kin{cond}.mean_speed(array_data{u}.stimData{cond});
            end
            
            keep_mask = trial_speeds > speed_bounds(1) & trial_speeds < speed_bounds(2);
            
            spike_trial_times = array_data{u}.spikeTrialTimes{cond}(keep_mask==1);
            array_data_trim{u}.spikeTrialTimes{cond} = spike_trial_times;
            array_data_trim{u}.binCounts{cond} = histcounts(spike_trial_times*1000,array_data_trim{u}.binEdges{cond});
            
            if(isfield(array_data{u},'num_stims')) % ugh
                array_data_trim{u}.num_stims(cond) = sum(array_data{u}.kin{cond}.mean_speed > speed_bounds(1) & array_data{u}.kin{cond}.mean_speed < speed_bounds(2));
                num_stims_cond = array_data_trim{u}.num_stims(cond);
            else
                array_data_trim{u}.numStims(cond) = sum(array_data{u}.kin{cond}.mean_speed > speed_bounds(1) & array_data{u}.kin{cond}.mean_speed < speed_bounds(2));
                num_stims_cond = array_data_trim{u}.numStims(cond);
            end
            
            if(normalizeBinCounts)
                array_data_trim{u}.binCounts{cond} = array_data_trim{u}.binCounts{cond}/num_stims_cond;
            end
            
            if(max(array_data_trim{u}.binCounts{cond}) > array_data_trim{u}.binMaxYLim)
               array_data_trim{u}.binMaxYLim = max(array_data_trim{u}.binCounts{cond})*1.1; 
            end
            
        end
    end




end

