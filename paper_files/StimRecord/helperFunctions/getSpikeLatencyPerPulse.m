function [output_data] = getSpikeLatencyPerPulse(array_data,pulse_width)

    if(~iscell(array_data))
        array_data = {array_data};
    end

    spike_times = cell(1,numel(array_data{1}.binCounts));
    
    % for each unit
    for unit_idx = 1:numel(array_data)
        % for each condition
        for condition = 1:numel(array_data{1}.binCounts)
            % initialize spike_times output if first unit
            if(unit_idx == 1)
                spike_times{condition} = cell(numel(array_data{unit_idx}.PULSE_TIMES{condition}{1}),1);
            end
            
            % for each train
            for train_idx = 1:array_data{unit_idx}.numStims(condition)
                % get spike latency relative to each pulse in the train
                window = [0,mode(diff(array_data{unit_idx}.PULSE_TIMES{condition}{train_idx}))]; % in s, remove the last 0.5ms
                spike_trial_times = array_data{unit_idx}.spikeTrialTimes{condition}(array_data{unit_idx}.stimData{condition} == train_idx);
                
                for pulse = 1:numel(array_data{unit_idx}.PULSE_TIMES{condition}{train_idx})
                    window_adj = window + array_data{unit_idx}.PULSE_TIMES{condition}{train_idx}(pulse); 

                    spike_pulse_mask = spike_trial_times > window_adj(1) & spike_trial_times <= window_adj(2);
                    spike_times{condition}{pulse} = [spike_times{condition}{pulse},spike_trial_times(spike_pulse_mask == 1) - window_adj(1) - pulse_width];

                end
            end
        end
    end


    output_data.spike_times = spike_times;


end