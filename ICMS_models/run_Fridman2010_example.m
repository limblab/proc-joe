%% runs Fridman model

    PARAMS.Ks = 1; % perceived intensity gain
    PARAMS.perception_tau = 0.5; % s
    PARAMS.incorporate_stim_time_constant = 1;
    PARAMS.run_controller = 0;
    CONSTANTS.dt = 0.4/1000; % s, to match the biphasic pulse length
    CONSTANTS.num_pulses_per_dt = round((0.4/1000)/(CONSTANTS.dt));
    CONSTANTS.t_max = 40; % s
    CONSTANTS.min_current = 10; % uA
%% compute perceived intensity against time for different amplitudes

    currents = [60]; % uA
    amp_test = [60];
    rates_found = [3.4];
    rates_use = rates_found; %interp1(amp_test,rate_found,currents,'linear','extrap');
    
    STIM_PARAMS.recovery_tau = 0; % seconds until stim_decay kicks in
    stim_decay_taus = 1./rates_use;
    STIM_PARAMS.pulse_width = 200; % us
    STIM_PARAMS.frequency = 300; % Hz
    STIM_PARAMS.duty_cycle = 0.5;
    STIM_PARAMS.duty_length = 0.2; % s
    move_time = 0.12;
    rt_lowest = 0.35 - move_time; % correct for time to produce movement
 
    figure(); hold on
    on_estimate = nan(numel(currents),1);
    off_estimate = nan(numel(currents),1);
    colors = [(0+(1:numel(currents))*1/numel(currents))',zeros(numel(currents),1),zeros(numel(currents),1)];
    for i_current = 1:numel(currents)
        STIM_PARAMS.current = currents(i_current);
        STIM_PARAMS.stim_decay_tau = stim_decay_taus(i_current);
        perception_data = computePerceivedIntensity(STIM_PARAMS, PARAMS, CONSTANTS);
        plot(perception_data.t_list,perception_data.Ps_list,'color',colors(i_current,:));
        if(i_current == 1)
%             threshold = perception_data.Ps_list(find(perception_data.t_list >= rt_lowest,1,'first'));
            plot([0,perception_data.t_list(end)],[threshold,threshold],'--','color',getColorFromList(1,1));
        end
        on_idx = find(perception_data.Ps_list > threshold,1,'first');
        on_estimate(i_current) = CONSTANTS.dt*on_idx;
        if(~isempty(find(perception_data.Ps_list(on_idx+ceil(0.02/CONSTANTS.dt):end) < threshold,1,'first')))
            off_estimate(i_current) = CONSTANTS.dt*(on_idx + ceil(0.02/CONSTANTS.dt) + find(perception_data.Ps_list(on_idx+ceil(0.02/CONSTANTS.dt):end) < threshold,1,'first'));
        end
    end
    xlabel('Time (s)')
    ylabel('Perceived intensity (au)')
    formatForLee(gcf)
    set(gca,'fontsize',14)
    
    figure(); 
    plot(currents,on_estimate+move_time)
    figure();
    plot(currents(~isnan(off_estimate)),off_estimate(~isnan(off_estimate)))

%% compute perceived intensity against time for different frequencies

    STIM_PARAMS.current = 60; % uA
    frequencies = [20,100,300]; % Hz
    stim_taus = 1./[0.01,0.19,3.6];
    STIM_PARAMS.tau_offset = 0.0; % seconds until stim_decay kicks in
    STIM_PARAMS.pulse_width = 200; % us
    STIM_PARAMS.duty_cycle = [0.5]; % fraction
    STIM_PARAMS.duty_length = [0.2]; % s
    
    rt_lowest = 0.35 - 0.12; % correct for time to produce movement
    figure(); hold on
    on_estimate = zeros(numel(frequencies),1);
    colors = [(0+(1:numel(frequencies))*1/numel(frequencies))',zeros(numel(frequencies),1),zeros(numel(frequencies),1)];
    for i_freq = 1:numel(frequencies)
        % use duty cycle
        STIM_PARAMS.stim_decay_tau = stim_taus(i_freq);
        STIM_PARAMS.frequency = frequencies(i_freq);
        perception_data = computePerceivedIntensity(STIM_PARAMS, PARAMS, CONSTANTS);
        plot(perception_data.t_list,perception_data.Ps_list,'color',colors(i_freq,:));
        if(i_freq == 1)
            threshold = perception_data.Ps_list(find(perception_data.t_list >= rt_lowest,1,'first'));
%             plot([0,perception_data.t_list(end)],[threshold,threshold],'--','color',getColorFromList(1,1));
        end
        on_estimate(i_freq) = CONSTANTS.dt*find(perception_data.Ps_list > threshold,1,'first');
    end
    figure(); 
    plot(frequencies,on_estimate)
    
%% compute perceived intensity against time for different duty cycles, match charge injected (frequency)

    STIM_PARAMS.current = 40; % uA
    stim_taus = 1./[0.15,0.22,0.5,1.6];
%     taus = [1,1,1,1]*(1/0.15);
    STIM_PARAMS.pulse_width = 200; % us
    base_frequency = 180; % Hz
    duty_cycles = [0.33,0.5,0.67,1];
    duty_lengths = [0.3,0.2,0.15,1]; % s
    
    figure(); hold on
    colors = [(0+(1:numel(duty_cycles))*1/numel(duty_cycles))',zeros(numel(duty_cycles),1),zeros(numel(duty_cycles),1)];
    for i_cycle = 1:numel(duty_cycles)
        % use duty cycle
        STIM_PARAMS.tau = stim_taus(i_cycle);
        STIM_PARAMS.duty_cycle = duty_cycles(i_cycle);
        STIM_PARAMS.duty_length = duty_lengths(i_cycle);
        STIM_PARAMS.frequency = base_frequency;
        perception_data = computePerceivedIntensity(STIM_PARAMS, PARAMS, CONSTANTS);
        plot(perception_data.t_list,perception_data.Ps_list,'color',colors(i_cycle,:));
        
        % use control frequency and continuous stimulation
%         STIM_PARAMS.tau = taus(i_cycle);
%         STIM_PARAMS.duty_cycle = 1;
%         STIM_PARAMS.duty_length = duty_lengths(i_cycle);
%         STIM_PARAMS.frequency = base_frequency*duty_cycles(i_cycle);
%         perception_data = computePerceivedIntensity(STIM_PARAMS, PARAMS, CONSTANTS);
%         plot(perception_data.t_list,perception_data.Ps_list,'-.','color',colors(i_cycle,:));
        
    end



%% runs Fridman model with controller?

    PARAMS.Ks = 1; % perceived intensity gain
    PARAMS.perception_tau = 0.5; % s
    PARAMS.incorporate_stim_time_constant = 1;
    PARAMS.run_controller = 1;
        
    CONSTANTS.controller_params.Kp = 150; 
    CONSTANTS.controller_params.Ki = 150; 
    CONSTANTS.controller_params.Kd = 5;
    CONSTANTS.desired_Ps = 1;
    
    CONSTANTS.dt = 0.4/1000; % s, to match the biphasic pulse length
    CONSTANTS.num_pulses_per_dt = round((0.4/1000)/(CONSTANTS.dt));
    
    CONSTANTS.t_max = 5; % s
    CONSTANTS.min_current = 0; % uA
%% compute perceived intensity against time for different amplitudes

    currents = [0]; % uA
    STIM_PARAMS.recovery_tau = 0; % seconds until stim_decay kicks in
    STIM_PARAMS.stim_decay_tau = 1./2;
    STIM_PARAMS.pulse_width = 200; % us
    STIM_PARAMS.frequency = 50; % Hz
    STIM_PARAMS.duty_cycle = 1;
    STIM_PARAMS.duty_length = 0.2; % s
    move_time = 0.12;
    rt_lowest = 0.35 - move_time; % correct for time to produce movement
 
    figure('Position',[494 527 1316 420]); hold on
    on_estimate = nan(numel(currents),1);
    off_estimate = nan(numel(currents),1);
    colors = [(0+(1:numel(currents))*1/numel(currents))',zeros(numel(currents),1),zeros(numel(currents),1)];
    for i_current = 1:numel(currents)
        STIM_PARAMS.current = currents(i_current);
        
        perception_data = computePerceivedIntensity(STIM_PARAMS, PARAMS, CONSTANTS);
        subplot(1,3,1); hold on
        plot(perception_data.t_list,perception_data.Ps_list,'color',colors(i_current,:)); 
        subplot(1,3,2); hold on
        plot(perception_data.t_list,perception_data.current_list,'color',colors(i_current,:)); 
        subplot(1,3,3); hold on
        plot(perception_data.integral_error); 
        plot(perception_data.derivative_error);
    end
