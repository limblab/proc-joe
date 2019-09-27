function [output_data] = getPeakIndependence(array_data, latency_struct, input_data)

    % given peaks and a window around the peaks, determine if firing in one
    % peak depends on the other peaks using probability rules
        
    output_data = [];
    spike_trials = {}; % list of trials with a spike for each amp and peak
    pair_data = {};
    % populate spike_trials using latencies in latency_struct and
    % array_data.spikeTrialTimes;
    % latency_struct.latencies is in ms, array_data.spikeTrialTimes is in s

    for a = 1:numel(latency_struct.latencies) % for each condition
        spike_trials{a} = cell(numel(latency_struct.latencies{a}),1);
        for p = 1:numel(spike_trials{a}) % for each peak
            spike_trials{a}{p} = array_data.stimData{a}(find(abs(array_data.spikeTrialTimes{a}*1000 - latency_struct.latencies{a}(p)) < input_data.PEAK_WIDTH));
        end
        
        % for peak pair in spike_trials{a}, see if prob spike depends on earlier peak
        % firing
        pair_data{a} = [];
        peak_pairs = [];
        if(numel(spike_trials{a}) > 1) % need at least 2 peaks
            peak_pairs = nchoosek(1:numel(spike_trials{a}),2); % n x 2 consisting of idx pairs to compare
            pair_data{a} = zeros(size(peak_pairs,1),2);
            peak_pairs = [1,2];
            for pair_idx = 1:size(peak_pairs,1)
                % compute prob of spike in 2nd peak given spike happened in 1st peak
                % and compute prob of spike in 2nd peak given no spike in 1st
                % peak
                num_stims_total = array_data.numStims(a);
                num_stims_1st_peak = unique(numel(spike_trials{a}{peak_pairs(pair_idx,1)}));
                
                % get Prob(2nd | 1st)
                pair_data{a}(pair_idx,1) = numel(intersect(spike_trials{a}{peak_pairs(pair_idx,1)},spike_trials{a}{peak_pairs(pair_idx,2)}))/num_stims_1st_peak;
                
                % get Prob(2nd)
                pair_data{a}(pair_idx,2) = (numel(spike_trials{a}{peak_pairs(pair_idx,2)})/num_stims_total);
                
                % if independent, these two things should be the same
                
                
                
            end
        end
    end

    output_data.pair_data = pair_data;
    output_data.spike_trials = spike_trials;

    
end