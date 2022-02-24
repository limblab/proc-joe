function [array_data] = getBaselineFiringRate(array_data,window)

    for i_unit = 1:numel(array_data)
        num_stims = sum(sum(array_data{i_unit}.numStims));
        num_spikes_window = 0;
        for i_amp = 1:numel(array_data{i_unit}.spikeTrialTimes)
            num_spikes_window = num_spikes_window + sum(array_data{i_unit}.spikeTrialTimes{i_amp} > window(1) & array_data{i_unit}.spikeTrialTimes{i_amp} < window(2));
        end
        array_data{i_unit}.baseline_fr = (num_spikes_window/abs(diff(window)*num_stims));
    end

end