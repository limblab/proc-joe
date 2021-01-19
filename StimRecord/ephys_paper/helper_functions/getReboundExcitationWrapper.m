function [output_data] = getReboundExcitationWrapper(array_data,input_data)

    filtered_PSTH = [];
    PSTH = [];
    amp = [];
    is_rebound = nan(numel(array_data),numel(input_data.cond_list));
    rebound_dur = nan(numel(array_data),numel(input_data.cond_list));
    rebound_mag = nan(numel(array_data),numel(input_data.cond_list));
    threshold = nan(numel(array_data),numel(input_data.cond_list));
    
    orig_post_window = input_data.post_window;
    orig_blank_time = input_data.blank_time;
    for i_unit = 1:numel(array_data)
        % for each condition:
        
        % sanitize array_data
        % santize array data

        if(~isfield(array_data{i_unit},'num_stims') && isfield(array_data{i_unit},'numStims'))
            array_data{i_unit}.num_stims = array_data{i_unit}.numStims;
        end
        
        for i_cond = 1:numel(input_data.cond_list)
            cond_idx = input_data.cond_list(i_cond);
            if(numel(array_data{i_unit}.num_stims) >= 8 && array_data{i_unit}.num_stims(cond_idx) > 0) % 
                % update post window and blank time to look at last pulse
                % in sequence
                input_data.post_window = orig_post_window + array_data{i_unit}.PULSE_TIMES{cond_idx}{4}(end);
                input_data.blank_time = orig_blank_time + [array_data{i_unit}.PULSE_TIMES{cond_idx}{4}(1),...
                    array_data{i_unit}.PULSE_TIMES{cond_idx}{4}(end)];
                
                input_data.cond_idx = cond_idx;
                temp_data = getReboundExcitation(array_data{i_unit},input_data);
                
                PSTH(i_unit,i_cond,:) = temp_data.PSTH;
                filtered_PSTH(i_unit,i_cond,:) = temp_data.filtered_PSTH;
                is_rebound(i_unit,i_cond) = temp_data.is_rebound;
                rebound_dur(i_unit,i_cond) = temp_data.rebound_dur;
                rebound_mag(i_unit,i_cond) = temp_data.rebound_mag;
                threshold(i_unit,i_cond) = temp_data.threshold;
                
                input_data.post_window = orig_post_window;
                input_data.blank_time = orig_blank_time;
            end
        end
    end
    
    output_data.filtered_PSTH = filtered_PSTH;
    output_data.is_rebound = is_rebound;
    output_data.rebound_dur = rebound_dur;
    output_data.rebound_mag = rebound_mag;
    output_data.PSTH = PSTH;
    output_data.threshold = threshold;
    output_data.amp = amp;


end