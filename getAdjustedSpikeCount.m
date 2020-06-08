function [trial_data] = getAdjustedSpikeCount(trial_data,params)
    % adjusts spike counts based on how many stimulations occur in each bin
    
    % default params
    IPI = 10; num_pulses = 10;
    stim_field = 'idx_stimTime';
    array_name = 'LeftS1';
    artifact_length = 1/1000;
    if nargin > 1, assignParams(who,params); end
    
    num_stims_per_bin = histcounts((0:(num_pulses-1))*IPI,'BinWidth',trial_data(1).bin_size*1000);
    
    for i_trial = 1:numel(trial_data)
        if(~isnan(trial_data(i_trial).idx_stimTime)) % we have a stimulation
            idx_stim = trial_data(i_trial).idx_stimTime;
            trial_data(i_trial).([array_name,'_spikes'])(idx_stim:idx_stim+numel(num_stims_per_bin)-1,:) = ...
                trial_data(i_trial).([array_name,'_spikes'])(idx_stim:idx_stim+numel(num_stims_per_bin)-1,:).*(trial_data(i_trial).bin_size./(trial_data(i_trial).bin_size-artifact_length*num_stims_per_bin))';
        end
    end
    



end