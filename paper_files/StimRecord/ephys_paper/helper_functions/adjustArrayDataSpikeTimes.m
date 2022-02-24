function [array_data] = adjustArrayDataSpikeTimes(array_data,wave_length)

    for i_unit = 1:numel(array_data)
        for i_cond = 1:numel(array_data{i_unit}.spikeTrialTimes)
            array_data{i_unit}.spikeTrialTimes{i_cond} = array_data{i_unit}.spikeTrialTimes{i_cond} - wave_length;
        end
    end


end