%% assumes analyzeStimData was ran and arrayData is in the workspace

%% get firing rate during stimulation for each condition
    baseline_window = [-150,-5];
    window = [0,110];
    
    baseline_firing_rate = zeros(size(arrayData{1}.spikeTrialTimes));
    stim_firing_rate = zeros(size(arrayData{1}.spikeTrialTimes));
    firing_rate_evoked = zeros(size(arrayData{1}.spikeTrialTimes));
    
    for u = 1:numel(arrayData)
        baseline_idx = [find(arrayData{u}.binEdges{1} > baseline_window(1),1,'first'), ...
            find(arrayData{u}.binEdges{1} > baseline_window(2),1,'first')];
        stim_idx = [find(arrayData{u}.binEdges{1} > window(1),1,'first'),...
            find(arrayData{u}.binEdges{1} > window(2),1,'first')];
        
        for chan = 1:size(arrayData{u}.spikeTrialTimes,1)
            for wave = 1:size(arrayData{u}.spikeTrialTimes,2)
                baseline_firing_rate(chan,wave) = sum(arrayData{u}.binCounts{chan,wave}(baseline_idx(1):baseline_idx(2)))/(diff(baseline_window))*1000;
                stim_firing_rate(chan,wave) = sum(arrayData{u}.binCounts{chan,wave}(stim_idx(1):stim_idx(2)))/(diff(window))*1000;
               
            end
        end
        
        firing_rate_evoked = stim_firing_rate - baseline_firing_rate;
 
    end




