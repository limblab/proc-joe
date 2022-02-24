function [output_data] = convertWaveMappings(wave_mapping, input_data)
    

    % get actual frequencies instead of normalized frequencies
    if(max(wave_mapping(:,2)) <=1)
        wave_mapping(:,2) = wave_mapping(:,2)*input_data.max_freq;
    end

    % setup put wave_mappings and input_data into correct format
    input_data.freqs = wave_mapping(:,2);
    input_data.chans = wave_mapping(:,1);
    % call makeStimPattern
    
    output_data = makeStimPattern(input_data);
    
    % make sure output_data.stim_pattern is a valid pattern
    if(~isempty(find(sum(output_data.stim_pattern,1) > 16)))
        warning('Not a valid stim pattern. OOOPS');
    end
    
end