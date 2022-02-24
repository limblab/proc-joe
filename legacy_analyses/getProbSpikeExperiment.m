function [ bin_counts_exp,num_neurons_exp, prob_spike ] = getProbSpikeExperiment( arrayData,spikesStruct,amps,bin_edges, do_mean )

    
    bin_counts_exp = nan(numel(amps),numel(bin_edges)-1,numel(arrayData));
    num_neurons_exp = zeros(numel(amps),1);
    prob_spike = zeros(numel(amps),numel(arrayData));
    for i_amp = 1:numel(amps)
        % bin experimental data
        for i_neuron = 1:numel(arrayData)
            amp_idx = find([arrayData{i_neuron}.STIM_PARAMETERS.amp1] == amps(i_amp),1,'first');
            amp_100uA_idx = find(spikesStruct{i_neuron}.amp == 100,1,'first');
            if(~isempty(amp_idx))
                baseline_count = spikesStruct{i_neuron}.mean_baseline_fr*mode(diff(bin_edges))/1000;
                bin_counts_exp(i_amp,:,i_neuron) =  histcounts(arrayData{i_neuron}.spikeTrialTimes{amp_idx} - 0.453/1000,bin_edges/1000)/arrayData{i_neuron}.numStims(amp_idx) - baseline_count;
                num_neurons_exp(i_amp) = num_neurons_exp(i_amp) + 1;
                
                prob_spike(i_amp,i_neuron) = histcounts(arrayData{i_neuron}.spikeTrialTimes{amp_idx} - 0.453/1000, [0,5]/1000)/arrayData{i_neuron}.numStims(amp_idx);
            end
        end
    end


    if(do_mean)
        bin_counts_exp = mean(bin_counts_exp,3,'omitnan');
    end

end

