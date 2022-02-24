function [ output_data ] = updatePerception( input_data )

    %unpack input_data
    i_t = input_data.i_t; % time idx
    Ps_list = input_data.Ps_list; % perception list
    current = input_data.current;
    stim_on = input_data.stim_on;
    min_current = input_data.min_current;
    Ks = input_data.Ks;
    incorporate_stim_time_constant = input_data.incorporate_stim_time_constant;
    dt = input_data.dt;
    recovery_tau = input_data.recovery_tau;
    stim_decay_tau = input_data.stim_decay_tau;
    stim_effect =  input_data.stim_effect;
    perception_tau = input_data.perception_tau;
    duty_cycle_on = input_data.duty_cycle_on;
    
    % update Ps
    if(i_t == 1) % incorporate previous intensity
        prev_ps = 0;
    else
        prev_ps = Ps_list(i_t-1);
    end
    % incorporate stim effect
    if(stim_on(i_t) && current > min_current)
        stim_ps = Ks*(current^(3/2) - min_current^(3/2));
    else
        stim_ps = 0;
    end
    % incorporate stim effect time constant
    if(incorporate_stim_time_constant)
        % update stim_effectiveness
            if(duty_cycle_on(i_t))
                stim_effect = max(0,stim_effect - stim_effect/stim_decay_tau*dt);
            else
                stim_effect = min(1,stim_effect + recovery_tau*dt);
            end
        stim_ps = stim_ps*stim_effect;
    end
    Ps_list(i_t) = prev_ps + (-1/perception_tau)*prev_ps*dt + stim_ps*dt;

    % package outputs
    output_data.Ps_list = Ps_list;
    output_data.stim_effect = stim_effect;
    
end

