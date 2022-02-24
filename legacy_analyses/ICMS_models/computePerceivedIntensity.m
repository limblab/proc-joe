function [ output_data ] = computePerceivedIntensity(STIM_PARAMS, PARAMS,CONSTANTS)
% this function computes the number of spikes for the given inputs

% STIM_PARAMS -- current, pulse_width, num_pulses, frequency, duty cycle,
% duty length

% PARAMS -- Ks, tau

% CONSTANTS -- dt, t_max, current_min


    initialize_data = initializePerceptionModel(STIM_PARAMS,PARAMS,CONSTANTS);
    update_data = [];
    
    % things that don't change as for each update
    update_data.stim_on = initialize_data.stim_on;
    update_data.duty_cycle_on = initialize_data.duty_cycle_on;
    update_data.min_current = CONSTANTS.min_current;
    update_data.Ks = PARAMS.Ks;
    update_data.incorporate_stim_time_constant = PARAMS.incorporate_stim_time_constant;
    update_data.dt = CONSTANTS.dt;
    update_data.recovery_tau = STIM_PARAMS.recovery_tau;
    update_data.perception_tau = PARAMS.perception_tau;
    update_data.stim_decay_tau = STIM_PARAMS.stim_decay_tau;
    update_data.current = STIM_PARAMS.current;
    
    % initialize params that will change in for loop
    update_data.Ps_list = initialize_data.Ps_list;
    update_data.stim_effect = 1;
    stim_effect_list = zeros(size(update_data.Ps_list));
    current_list = zeros(size(update_data.Ps_list));
    error = []; derivative_error = []; integral_error = [];
    % update perceived intensity at each time point
    for i_t = 1:numel(update_data.Ps_list)
        update_data.i_t = i_t; % time idx
        stim_effect_list(i_t) = update_data.stim_effect;
        current_list(i_t) = update_data.current;
        perception_data = updatePerception(update_data);
                
        % update some variables
        update_data.Ps_list = perception_data.Ps_list; % perception list        
        update_data.stim_effect =  perception_data.stim_effect;
        
        % update current using a PID controller
        if(PARAMS.run_controller && update_data.stim_on(i_t)) % we are stimulating
            controller_ps_list = update_data.Ps_list(CONSTANTS.num_pulses_per_dt:initialize_data.stim_diff:end);
            controller_dt = initialize_data.stim_diff*CONSTANTS.dt;
            controller_idx = ceil(i_t/initialize_data.stim_diff);
            [update_data.current,errors] = runPerceptionController(CONSTANTS.controller_params,CONSTANTS.desired_Ps,controller_ps_list,controller_idx,controller_dt);
            error(end+1,1) = errors.error;
            integral_error(end+1,1) = errors.integral_error;
            derivative_error(end+1,1) = errors.derivative_error;
        end

    end
    
    % format output_data
    output_data.Ps_list = update_data.Ps_list;
    output_data.t_list = initialize_data.t_list;
    output_data.stim_on = update_data.stim_on;
    output_data.stim_times = initialize_data.stim_times;
    output_data.stim_effect_list = stim_effect_list;
    output_data.current_list = current_list;
    output_data.error = error;
    output_data.integral_error = integral_error;
    output_data.derivative_error = derivative_error;
end

