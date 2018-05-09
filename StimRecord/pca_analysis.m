%% convert to trial data
params = [];
params.include_start = 1;
params.extra_time = [0.2, 0.2];
trial_data = parseFileByTrial(cds,params);

%% add stim data to td -- this might be wrong, also not using cds.stimOn which is stupid of me
for td = 1:numel(trial_data)
    trial_data(td).stimTime = cds.trials.stimTime(trial_data(td).trial_id) - cds.trials.startTime(trial_data(td).trial_id) + params.extra_time(1);
    trial_data(td).idx_stimTime = ceil(trial_data(td).stimTime*100);
    trial_data(td).stimCode = cds.trials.stimCode((trial_data(td).trial_id));
    nonStim_trials(td) = isnan(trial_data(td).stimTime);
end

%% get pca data for non-stim trials
params_PCA.signals = {'LeftS1_spikes'};
params_PCA.trial_idx = nonStim_trials;
params_PCA.do_plot = 0;
[trial_data,pca_info] = getPCA(trial_data, params_PCA);

%% look at pc space
trial_idx = 3;

test = trial_data(trial_idx).LeftS1_spikes*pca_info.w;

%% visData
params_vD.trials_to_plot = 1;
params_vD.signals = {'pos'};
params_vd.continuous_data = 'pos';
params_vD.ARRAY_spikes = 'LeftS1';
params_vD.pca_array = 'LeftS1';
params_vD.plot_pca = 1;
visData(trial_data,params_vD);


