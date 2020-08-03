function [output_data] = getInhibitionDurationDoubplePulseWrapper(array_data,input_data)

    % gets inhibition duration for the i_ampitions in array_data. Inhibition
    % duration is defined based on (Butovas 2003)
    
    % input_data contains
    %   PRE_WINDOW window used to compute spontaneous firing rate in ms
    %   POST_WINDOW window used to compute inhib duration (if exists)
    %   BLANK_TIME -- time after stimulation onset to blank to remove
    %       excitatory period due to stimulation
      
    filtered_PSTH = [];
    PSTH = [];
    amp = [];
    threshold = [];
    is_inhib = zeros(numel(array_data),numel(input_data.cond_list),1);
    inhib_dur = nan(numel(array_data),numel(input_data.cond_list),1);
    threshold = nan(numel(array_data),numel(input_data.cond_list),1);
    
    orig_post_window = input_data.post_window;
    orig_blank_time = input_data.blank_time;
    orig_time_start = input_data.max_time_start;
    
    for i_unit = 1:numel(array_data)
        % for each i_ampition:
        for i_cond = 1:numel(input_data.cond_list)
            cond_idx = input_data.cond_list(i_cond);
            if(numel(array_data{i_unit}.num_stims) == 8 && array_data{i_unit}.num_stims(cond_idx) > 0) % 
                
                % update post window and blank time to look at last pulse
                % in sequence
%                 if(i_cond ~= find(input_data.cond_list == 6,1,'first'))
                    input_data.post_window = input_data.post_window + array_data{i_unit}.PULSE_TIMES{cond_idx}{1}(end);
                    input_data.blank_time = input_data.blank_time + array_data{i_unit}.PULSE_TIMES{cond_idx}{1}(end);
%                 end
                
                temp_inhib_data = getInhibitionDuration(array_data{i_unit},cond_idx,input_data);
                
                PSTH(i_unit,i_cond,:) = temp_inhib_data.PSTH;
                filtered_PSTH(i_unit,i_cond,:) = temp_inhib_data.filtered_PSTH;
                is_inhib(i_unit,i_cond) = temp_inhib_data.is_inhib;
                inhib_dur(i_unit,i_cond) = temp_inhib_data.inhib_dur;
                threshold(i_unit,i_cond) = temp_inhib_data.threshold;
                
                input_data.post_window = orig_post_window;
                input_data.blank_time = orig_blank_time;
                input_data.max_time_start = orig_time_start;
            end
        end
    end
    
    output_data.filtered_PSTH = filtered_PSTH;
    output_data.is_inhib = is_inhib;
    output_data.inhib_dur = inhib_dur;
    output_data.PSTH = PSTH;
    output_data.threshold = threshold;
    output_data.amp = amp;
end