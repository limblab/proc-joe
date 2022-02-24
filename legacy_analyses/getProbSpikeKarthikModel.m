function [ bin_counts, num_spikes_post_stim ] = getProbSpikeKarthikModel( neuron_struct, bin_edges, dist_range, num_stims,amps, wave_length, do_mean )
    % bin model data at a specified bin edge
    % returns either an amp x bin edge matrix if do_mean == 1
    % or an amp x bin edge x neuron matrix if do_mean == 0
    
    bin_counts = nan(numel(amps),numel(bin_edges)-1,numel(neuron_struct));
    num_neurons = zeros(numel(amps),1);
    num_spikes_post_stim = nan(numel(amps),numel(neuron_struct));
    
    for i_amp = 1:numel(amps)
        % bin model data
        for i_neuron = 1:numel(neuron_struct)
            amp_idx = find(neuron_struct{i_neuron}.amps == amps(i_amp));
            if(~isempty(amp_idx) && neuron_struct{i_neuron}.dist > dist_range(1) && neuron_struct{i_neuron}.dist < dist_range(2))
                baseline_count = mean(neuron_struct{i_neuron}.baseline_fr)*mode(diff(bin_edges))/1000;
                bin_counts(i_amp,:,i_neuron) = histcounts(neuron_struct{i_neuron}.spike_times{amp_idx} - wave_length,bin_edges)/num_stims - baseline_count;
                num_neurons(i_amp) = num_neurons(i_amp) + 1;
                
                num_spikes_post_stim(i_amp,i_neuron) = histcounts(neuron_struct{i_neuron}.spike_times{amp_idx} - wave_length,[0,5]);
            end
        end        
        
    end
    
    if(do_mean == 1)
        bin_counts = sum(bin_counts,3,'omitnan')./num_neurons;
    end

end

