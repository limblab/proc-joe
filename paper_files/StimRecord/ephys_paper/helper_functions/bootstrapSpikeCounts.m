function [ spike_count_all ] = bootstrapSpikeCounts( input_data )

    spike_count_all = [];
    for i = 1:input_data.n
        spike_times = datasample(input_data.base_spike_times,numel(input_data.base_spike_times),'Replace',true);
        spike_count = histcounts(spike_times,input_data.bin_edges);
        
        spike_count_all = [spike_count_all,spike_count];
    end
    
    spike_count_all = spike_count_all/input_data.n_trials;
end

