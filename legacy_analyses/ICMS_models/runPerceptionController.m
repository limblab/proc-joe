function [stim_param, errors] = runPerceptionController(controller, desired_Ps, Ps_list, i_Ps, dt)

    % controller params has Kp, Kd, Ki (constants for PID controller)
    % controller params

    Kp = controller.Kp; Kd = controller.Kd; Ki = controller.Ki;
    error = desired_Ps - Ps_list(i_Ps);
    integral_error = sum((desired_Ps-Ps_list(1:i_Ps))*dt);
    if(i_Ps > 2)
        derivative_error = -1*(Ps_list(i_Ps) - Ps_list(i_Ps-2))/(2*dt);
    else
        derivative_error = 0;
    end
    
    stim_param = max(0,Kp*error + Ki*integral_error + Kd*derivative_error);

    errors.error = error;
    errors.integral_error = integral_error;
    errors.derivative_error = derivative_error;
end