function [spikes_in_window] = getSpikesInWindow(array_data,condition_compute,window,output_fr)

    if(isempty(condition_compute))
        condition_compute = 1:1:numel(array_data.binCounts);
    end

    if(~isfield(array_data,'num_stims'))
        array_data.num_stims = array_data.numStims;
    end
    
    spikes_in_window = zeros(1,numel(condition_compute));
    spikes_in_window_idx = 1;
    for condition = condition_compute
        spikes_in_window(spikes_in_window_idx) = sum(array_data.spikeTrialTimes{condition}*1000 > window(1) & ...
            array_data.spikeTrialTimes{condition}*1000 < window(2))/array_data.num_stims(condition);
        spikes_in_window_idx = spikes_in_window_idx + 1;
    end
    
        
    % convert to fr
    if(output_fr)
        spikes_in_window = spikes_in_window/(diff(window))*1000; % convert to Hz
    end



end