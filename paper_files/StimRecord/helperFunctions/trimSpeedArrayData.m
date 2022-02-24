function [array_data_trim] = trimSpeedArrayData(array_data,mean_speed_cutoff)

    % returns array_data_trim with stim 'trials' where the mean speed is
    % between the cutoffs in mean_speed_cutoff
    
    
    array_data_trim = array_data;
    for u = 1:numel(array_data)
        for cond = 1:numel(array_data{u}.binCounts)
            trial_keep_mask = array_data{u}.trial_num{cond} < numel(array_data{u}.kin{cond}.mean_speed);
            trial_speeds = array_data{u}.kin{cond}.mean_speed(array_data{u}.trial_num{cond}(trial_keep_mask));
            keep_mask = trial_speeds > mean_speed_cutoff(1) & trial_speeds < mean_speed_cutoff(2);
            trials_keep = unique(array_data{u}.trial_num{cond}(keep_mask));


            spike_trial_times = array_data{u}.spikeTrialTimes{cond}(keep_mask==1);
            array_data_trim{u}.binCounts{cond} = histcounts(spike_trial_times*1000,array_data_trim{u}.binEdges{cond});
            array_data_trim{u}.num_stims(cond) = sum(array_data{u}.kin{cond}.mean_speed > mean_speed_cutoff(1) & array_data{u}.kin{cond}.mean_speed < mean_speed_cutoff(2));
            array_data_trim{u}.kin{cond}.speed = array_data{u}.kin{cond}.speed(trials_keep,:);
            array_data_trim{u}.kin{cond}.mean_speed = array_data{u}.kin{cond}.mean_speed(trials_keep,:);
        end
    end


end