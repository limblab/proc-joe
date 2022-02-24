function [output_data] = getDecayShortTrains(array_data,input_data)

    slope = nan(numel(array_data),numel(input_data.cond_list));
    
    array_data = getResponseToEachPulse(array_data,input_data);
    
    for i_unit = 1:numel(array_data)
        for i_cond = 1:numel(input_data.cond_list)
            cond_idx = input_data.cond_list(i_cond);
            if(numel(array_data{i_unit}.num_stims) >= 8 && array_data{i_unit}.num_stims(cond_idx) > 0) % 
                % update post window and blank time to look at last pulse
                % in sequence
                
                input_data.cond_idx = cond_idx;
                temp_data = getDecay(array_data{i_unit},input_data);
                
                slope(i_unit,i_cond) = temp_data.slope; % FR/s
                beg_vs_end(i_unit,i_cond) = temp_data.beg_vs_end; % FR
            end
        end
        
        
    end

    % package outputs
    output_data.slope = slope;
    output_data.beg_vs_end = beg_vs_end;
end


function [output_data] = getDecay(data,input_data)
    % computes two metrics:
    % slope of a line fit through mean response to each pulse
    % difference between mean response of first and last pulse

    m = nan; beg_vs_end = nan;
    
    % fit line
    mdl = fit(data.PULSE_TIMES{input_data.cond_idx}{1}',data.mean_resp_each_pulse{input_data.cond_idx}','m*x+b');
    m = mdl.m;
    
    % first vs. last;
    first = data.mean_resp_each_pulse{input_data.cond_idx}(1);
    last = data.mean_resp_each_pulse{input_data.cond_idx}(end);
    beg_vs_end = (first-last)/first;
    
    
    output_data.slope = m; % m is FR per s
    output_data.beg_vs_end = beg_vs_end;
end




