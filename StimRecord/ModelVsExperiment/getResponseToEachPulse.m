function [array_data] = getResponseToEachPulse(array_data, input_data)

    
    for i_unit = 1:numel(array_data)
        array_data{i_unit}.response_each_pulse = cell(size(array_data{i_unit}.stimData));
        
        for i_cond = 1:numel(array_data{i_unit}.response_each_pulse)
            if(array_data{i_unit}.numStims(i_cond) > 0)
                % initialize matrices
                array_data{i_unit}.response_each_pulse{i_cond} = zeros(array_data{i_unit}.numStims(i_cond),numel(array_data{i_unit}.PULSE_TIMES{i_cond}{1}));
                array_data{i_unit}.mean_resp_each_pulse{i_cond} = zeros(1,numel(array_data{i_unit}.PULSE_TIMES{i_cond}));

                % grab spike data
                spike_times = array_data{i_unit}.spikeTrialTimes{i_cond};
                trial_num = array_data{i_unit}.stimData{i_cond};
                num_pulses = numel(array_data{i_unit}.PULSE_TIMES{i_cond}{1});

                % get FR to each pulse
                for i_trial = 1:array_data{i_unit}.numStims(i_cond)
                    % get response for each pulse
                    for i_pulse = 1:num_pulses
                        offset = array_data{i_unit}.PULSE_TIMES{i_cond}{1}(i_pulse);
                        spike_mask = trial_num == i_trial & spike_times > input_data.spike_window(1)/1000 + offset & spike_times < input_data.spike_window(2)/1000 + offset;
                        array_data{i_unit}.response_each_pulse{i_cond}(i_trial,i_pulse) = sum(spike_mask)/diff(input_data.spike_window/1000);
                    end
                end

                % store mean FR for each pulse
                array_data{i_unit}.mean_resp_each_pulse{i_cond} = mean(array_data{i_unit}.response_each_pulse{i_cond},1);
                
            end
        end
   
    end



end