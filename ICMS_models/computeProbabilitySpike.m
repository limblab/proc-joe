function [prob,threshold_current] = computeProbabilitySpike(STIM_PARAMS, PARAMS, CONSTANTS)

% this function computes the probability of a spike for a given stim pulse

% STIM_PARAMS -- current, pulse_width, num_pulses, frequency

% PARAMS -- gain (for this task), spread, refractory_time_constant, conductance,...
% chronaxie, rheobase, exp_time_constant

% CONSTANTS -- t_abs, dt, distances
    
    % compute Iz
    threshold_current = computeThresholdCurrent(STIM_PARAMS,PARAMS,CONSTANTS)';
    Iz = (PARAMS.gain*STIM_PARAMS.current./(CONSTANTS.distances.^2) - threshold_current)./(PARAMS.spread*threshold_current);
    
    % compute probability
    prob = normcdf(Iz,0,1);

end

function [threshold_current] = computeThresholdCurrent(STIM_PARAMS,PARAMS,CONSTANTS)
    % computes the threshold current based on spike history

    threshold_current = PARAMS.rheobase*(1+PARAMS.chronaxie/STIM_PARAMS.pulse_width)*...
            (1+PARAMS.conductance*exp(-(CONSTANTS.dt-CONSTANTS.t_abs)/PARAMS.refractory_time_constant));
   
    threshold_current(CONSTANTS.dt < CONSTANTS.t_abs) = 10000;
end