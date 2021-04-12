function [output_data] = getBaselineActivityBlocked(array_data, input_data)
    output_data = [];

    for i_unit = 1:numel(array_data)
        for i_cond = 1:numel(input_data.cond_list)
            cond_idx = input_data.cond_list(i_cond);
            if(numel(array_data{i_unit}.num_stims) >= 8 && array_data{i_unit}.num_stims(cond_idx) > 0) % 
                blank_times = array_data{i_unit}.PULSE_TIMES{cond_idx}{4}(:) + [0,input_data.blank_time];
                num_baseline_spikes_expected = (input_data.max_train_time - input_data.blank_time*size(blank_times,2))*array_data{i_unit}.baseline_fr;
                % get spike during train
                spike_mask = array_data{i_unit}.spikeTrialTimes{cond_idx} > 0 & array_data{i_unit}.spikeTrialTimes{cond_idx} < input_data.max_train_time;
                % remove those within blanking windows
                for i_blank = 1:size(blank_times,2)
                    spike_mask = spike_mask & ~(array_data{i_unit}.spikeTrialTimes{cond_idx} >= blank_times(i_blank,1) & array_data{i_unit}.spikeTrialTimes{cond_idx} < blank_times(i_blank,2));
                end
                
                spike_ts = array_data{i_unit}.spikeTrialTimes{cond_idx}(spike_mask);
                num_baseline_spikes = numel(spike_ts)/array_data{i_unit}.num_stims(cond_idx);
                
                
                
                
            end
        end
    end


end