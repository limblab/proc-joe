function [ output_data ] = computeNumberSpikes(STIM_PARAMS, PARAMS, CONSTANTS)
% this function computes the number of spikes for the given inputs

% STIM_PARAMS -- current, pulse_width, num_pulses, frequency

% PARAMS -- gain (for this task), spread, refractory_time_constant, conductance,...
% chronaxie, rheobase, exp_time_constant

% CONSTANTS -- t_abs, dt, distances
    CONSTANTS.dt = 10000+zeros(numel(CONSTANTS.distances), 1);
    spike_history = zeros(numel(CONSTANTS.distances), STIM_PARAMS.num_pulses);
    prob_spike = zeros(numel(CONSTANTS.distances), STIM_PARAMS.num_pulses);
    threshold_current = zeros(numel(CONSTANTS.distances), STIM_PARAMS.num_pulses);
    
    for n = 1:STIM_PARAMS.num_pulses
            [prob_spike(:,n),threshold_current(:,n)] = computeProbabilitySpike(STIM_PARAMS, PARAMS, CONSTANTS);
            spike_history(:,n) = rand(size(prob_spike(:,n))) < prob_spike(:,n);
%             spike_history(:,n) = 0.5 > prob_spike(:,n);
            CONSTANTS.dt(spike_history(:,n) == 1) = 1000/STIM_PARAMS.frequency;
            CONSTANTS.dt(spike_history(:,n) == 0) = CONSTANTS.dt(spike_history(:,n) == 0) + 1000/STIM_PARAMS.frequency;
    end
    
    w_n = exp(-1000*(1:STIM_PARAMS.num_pulses)/STIM_PARAMS.frequency/PARAMS.exp_time_constant);
    num_spikes_total = sum(4*pi.*CONSTANTS.distances.^2.*sum(w_n.*prob_spike,2)');
    var_spikes_total = sum(4*pi.*CONSTANTS.distances.^2.*sum(w_n.*prob_spike.*(1-prob_spike),2)');
    
    
    % format output_data
    output_data.num_spikes_total = num_spikes_total;
    output_data.prob_spike = prob_spike;
    output_data.spike_history = spike_history;
    output_data.threshold_current = threshold_current;
    output_data.var_spikes_total = var_spikes_total;
end

