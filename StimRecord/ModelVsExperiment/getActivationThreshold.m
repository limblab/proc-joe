function [output_data] = getActivationThreshold(array_data, input_data)
% gets activation threshold (first amplitude that evokes a spike >= 50% of
% the time) for each entry in array_data


    thresholds = nan(numel(array_data),1);
    is_responsive = zeros(numel(array_data),1);
    
    for i_unit = 1:numel(array_data)
        % count # stims with spikes post stim in input_data.window
        percent_responsive = nan(numel(array_data{i_unit}.stimData),1);
        for i_amp = 1:numel(array_data{i_unit}.numStims)
            num_responsive = 0;
            total_stims = 0;
            for i_stim = 1:array_data{i_unit}.numStims(i_amp)
                % if we are either not removing intrinsic, the data is not
                % from the model, or intrinsic spikes are not near the
                % stim, continue.
                if(input_data.remove_intrinsic == 0 || input_data.is_model == 0 || array_data{i_unit}.intrinsic_near_stim{i_amp}(i_stim) == 0) 
                    spike_times = array_data{i_unit}.spikeTrialTimes{i_amp}(array_data{i_unit}.stimData{i_amp} == i_stim);
                    num_responsive = num_responsive + any(spike_times > input_data.spike_window(1) & spike_times <= input_data.spike_window(2));
                    total_stims = total_stims + 1;                    
                end
            end
            percent_responsive(i_amp) = num_responsive/total_stims;
        end
        
        unit_thresh = find(percent_responsive >= 0.5,1,'first');
        if(~isempty(unit_thresh))
            thresholds(i_unit) = array_data{i_unit}.STIM_PARAMETERS(unit_thresh).amp1;
            is_responsive(i_unit) = 1;
        end
    end

    % package outputs

    output_data.thresholds = thresholds;
    output_data.is_responsive = is_responsive;

end