function [output_data] = getBinnedSpikeTimes(array_data,input_data)
% input_data contains
    % window : s, window to bin spikes in
    % bin_size : s, size of bins
    % amp_list : list of amplitudes to be analyzed
    
    
    % make bins and bin_data
    bin_edges = input_data.window(1):input_data.bin_size:input_data.window(2);
    
    bin_data = [];
    for i_amp = 1:numel(input_data.amp_list)
        % meta info
        bin_data(i_amp).amp = input_data.amp_list(i_amp);
        bin_data(i_amp).diam_list = [];
        bin_data(i_amp).cell_id_list = [];
        bin_data(i_amp).bin_edges = bin_edges;
        bin_data(i_amp).bin_centers = bin_edges(1:end-1) + mode(diff(bin_edges))/2;
        
        % populate bin_data
        bin_data(i_amp).binned_data = []; % populate later
        for i_unit = 1:numel(array_data)
            amp_idx = find([array_data{i_unit}.STIM_PARAMETERS.amp1] == input_data.amp_list(i_amp),1,'first');
            if(~isempty(amp_idx)) % amp exists
                spike_mask = array_data{i_unit}.spikeTrialTimes{amp_idx} > input_data.window(1) & ...
                    array_data{i_unit}.spikeTrialTimes{amp_idx} <= input_data.window(2);
                
                sub_val = 0;
                if(input_data.sub_baseline)
                    sub_val = array_data{i_unit}.baseline_fr*input_data.bin_size;
                end
                
                bin_data(amp_idx).binned_data(end+1,:) = histcounts(array_data{i_unit}.spikeTrialTimes{1,amp_idx}(spike_mask),bin_edges)/array_data{i_unit}.numStims(amp_idx) - sub_val; % in # spikes per bin per stim
                
                if(input_data.is_model)
                    bin_data(amp_idx).diam_list(end+1,:) = array_data{i_unit}.diam;
                    bin_data(amp_idx).cell_id_list(end+1,:) = array_data{i_unit}.cell_id;
                end
            end
        end
    end
    

    output_data = bin_data;



end