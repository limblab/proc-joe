function [output_data] = getInhibitionDurationAmpWrapper(array_data,input_data)

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
    is_inhib = zeros(numel(array_data),numel(input_data.amp_list),1);
    inhib_dur = nan(numel(array_data),numel(input_data.amp_list),1);
        
    for i_unit = 1:numel(array_data)
        
        % for each i_ampition:
        for i_amp = 1:numel(input_data.amp_list)
        % low pass filter PSTH with a gaussian kernal, length 50 bins (1 ms bin
        % width)
            amp_idx = find([array_data{i_unit}.STIM_PARAMETERS.amp1] == input_data.amp_list(i_amp),1,'first');
            if(~isempty(amp_idx))
                
                temp_inhib_data = getInhibitionDuration(array_data{i_unit},amp_idx,input_data);
                PSTH(i_unit,i_amp,:) = temp_inhib_data.PSTH;
                filtered_PSTH(i_unit,i_amp,:) = temp_inhib_data.filtered_PSTH;
                is_inhib(i_unit,i_amp) = temp_inhib_data.is_inhib;
                inhib_dur(i_unit,i_amp) = temp_inhib_data.inhib_dur;
                
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



function data_smooth = gaussianKernel(data,kernel_length)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nbr_samples = numel(data);
    % apply smoothing to data

    % kernel half length is 3·SD out
    kernel_hl = ceil( 3 * kernel_length  );
    % create the kernel --it will have length 2*kernel_hl+1
    kernel = normpdf( -kernel_hl : ...
        1 : kernel_hl, ...
        0, kernel_length );
    
    % compute normalization factor --this factor depends on the number of taps
    % actually used
    nm = conv(kernel,ones(1,nbr_samples));

    % do the smoothing
        aux_smoothed_FR     = conv(kernel,data) ./ nm;
        % cut off the edges so that the result of conv is same length as the
        % original data
        data_smooth    = aux_smoothed_FR(kernel_hl+1:end-kernel_hl);
end
