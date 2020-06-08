function [ output_data ] = initializePerceptionModel( STIM_PARAMS, PARAMS, CONSTANTS )

    t = 0; % current time
    Ps_list = zeros(ceil(CONSTANTS.t_max/(CONSTANTS.dt)),1); % list of perceived intensities;
    t_list = 0:CONSTANTS.dt:CONSTANTS.t_max;
    t_list = t_list(1:numel(Ps_list)); % make sure these two matrices are the same size
    % get list of when stimulations occur
    stim_times = ((1:1:STIM_PARAMS.frequency*CONSTANTS.t_max)-1)*1/STIM_PARAMS.frequency;
    
    % make sure each pulse has the same number of time points
    num_time_points = CONSTANTS.num_pulses_per_dt;
    stim_on_idx = [1,find(diff(any(t_list-stim_times' <= STIM_PARAMS.pulse_width*2*1E-6 & t_list-stim_times' >=0))>.5)]; % first pulse starts immediately
    stim_diff = mode(diff(stim_on_idx));
    stim_on_idx = 1:stim_diff:numel(Ps_list);
    
    % adjust stim_on_idx based on the num of points for each pulse
    stim_on_idx = repmat(stim_on_idx,num_time_points,1); % for each time point
    stim_on_idx = stim_on_idx + [0:1:(num_time_points-1)]'; % offset
    stim_on_idx = reshape(stim_on_idx,numel(stim_on_idx),1);
    stim_on = zeros(size(Ps_list));
    stim_on(stim_on_idx) = 1;
       
    % adjust stim_on_idx based on duty cycle
    duty_length_idx = floor(STIM_PARAMS.duty_length/CONSTANTS.dt);
    duty_cycle_on = [ones(floor(duty_length_idx*STIM_PARAMS.duty_cycle),1);zeros(floor(duty_length_idx*(1-STIM_PARAMS.duty_cycle)),1)];
    duty_cycle_on = repmat(duty_cycle_on,ceil(numel(stim_on)/numel(duty_cycle_on)),1);
    stim_on = stim_on & duty_cycle_on(1:numel(stim_on));
    
    % package outputs
    output_data.Ps_list = Ps_list;
    output_data.t_list = t_list;
    output_data.stim_on = stim_on;
    output_data.stim_times = stim_times;
    output_data.duty_cycle_on = duty_cycle_on;
    output_data.stim_diff = stim_diff;
end

