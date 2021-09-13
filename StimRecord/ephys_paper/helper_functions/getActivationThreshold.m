function [output_data] = getActivationThreshold(array_data, input_data)
% gets activation threshold (first amplitude that evokes a spike >= 50% of
% the time) for each entry in array_data


    thresholds = nan(numel(array_data),1);
    is_responsive = zeros(numel(array_data),1);
    
    percent_responsive = nan(numel(array_data),numel(input_data.amp_list));
    percent_responsive_null = nan(numel(array_data),1);
    avg_num_spikes = nan(size(percent_responsive));
    avg_num_spikes_null = nan(size(percent_responsive_null));
    
    stat = nan(numel(array_data),numel(input_data.amp_list));
    
    for i_unit = 1:numel(array_data)
        % count # stims with spikes post stim in input_data.window
        for i_amp = 1:numel(array_data{i_unit}.numStims)
            amp_idx = find(input_data.amp_list == array_data{i_unit}.STIM_PARAMETERS(i_amp).amp1);
            if(~isempty(amp_idx))
                num_responsive = 0;
                total_stims = 0;
                num_spikes = 0;
                for i_stim = 1:array_data{i_unit}.numStims(i_amp)
                    % if we are either not removing intrinsic, the data is not
                    % from the model, or intrinsic spikes are not near the
                    % stim, continue.
                    if(input_data.remove_intrinsic == 0 || input_data.is_model == 0 || array_data{i_unit}.intrinsic_near_stim{i_amp}(i_stim) == 0) 
                        spike_times = array_data{i_unit}.spikeTrialTimes{i_amp}(array_data{i_unit}.stimData{i_amp} == i_stim);
                        num_responsive = num_responsive + any(spike_times > input_data.spike_window(1) & spike_times <= input_data.spike_window(2));
                        num_spikes = num_spikes + sum(spike_times > input_data.spike_window(1) & spike_times <= input_data.spike_window(2));
                        total_stims = total_stims + 1;                    
                    end
                end
                baseline_val = 0;
                if(input_data.use_stats || (input_data.sub_baseline && (input_data.is_model == 0 || input_data.remove_intrinsic == 0)))
                    baseline_val = array_data{i_unit}.baseline_fr*(diff(input_data.spike_window)); 
                end
                
%                 stat(i_unit,amp_idx) = poisscdf(num_spikes, baseline_val*total_stims,'upper');
                p_base = baseline_val; p_stim = num_responsive/total_stims; n_base = total_stims; n_stim = total_stims;
                stat(i_unit,amp_idx) = pearsonChiSquareBinomial(p_base,p_stim,n_base,n_stim,'lower');
                percent_responsive(i_unit,amp_idx) = p_stim;
                percent_responsive_null(i_unit) = p_base;
                avg_num_spikes(i_unit,amp_idx) = num_spikes/total_stims - baseline_val;
                avg_num_spikes_null(i_unit) = baseline_val;
            end
        end
        
        if(input_data.use_stats)
            unit_thresh = find(stat(i_unit,:) < input_data.threshold,1,'first');
        else
            unit_thresh = find(percent_responsive(i_unit,:) >= input_data.threshold,1,'first');
        end
        
        if(~isempty(unit_thresh))
            thresholds(i_unit) = input_data.amp_list(unit_thresh);
            is_responsive(i_unit) = 1;
        end
    end
    
    % package outputs

    output_data.thresholds = thresholds;
    output_data.is_responsive = is_responsive;
    output_data.percent_responsive = percent_responsive;
    output_data.percent_responsive_null = percent_responsive_null;
    output_data.poiss_stat = stat;
    output_data.avg_num_spikes_null = avg_num_spikes_null;
    output_data.avg_num_spikes = avg_num_spikes;
    
end