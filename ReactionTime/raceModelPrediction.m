%% define a population of exponential decreasing curves 
% RT = a*exp(-b*I)+c
    a_bounds = [0.1,0.25];
    b_bounds = [0.01,0.05];
    c_bounds = [0.15,0.17];
    std_bounds = [0.005,0.05];
    num_chans = 1000;

    a_all = rand(num_chans,1)*(diff(a_bounds)) + a_bounds(1);
    b_all = rand(num_chans,1)*(diff(b_bounds)) + b_bounds(1);
    c_all = rand(num_chans,1)*(diff(c_bounds)) + c_bounds(1);
    std_all = rand(num_chans,1)*diff(std_bounds) + std_bounds(1);

%% plot example RT vs amp curves
    figure();
    I_data = 1:1:100;
    RT_all = a_all.*exp(-b_all.*I_data)+c_all;

    plot(I_data,RT_all')
%% run my linear summation experiment
% sample N electrodes, compute RT based on the race model (sample each
% distribution, pick fastest). Store. Do this for different charges on each
% electrode and for different number of electrodes

    num_elecs = [1:24];
    total_charge = [240:120:600];
    num_runs_per_condition = 10000;

    RT_out = zeros(numel(num_elecs),numel(total_charge),num_runs_per_condition);
    for e = 1:numel(num_elecs)
        for c = 1:numel(total_charge)
            % sample num_elecs(e)
            chan_idx = [];
            for n = 1:num_runs_per_condition
                chan_idx(n,:) = randperm(num_chans,num_elecs(e));
            end
            % sample RT from normal distribution based on channel parameters
            a_sample = a_all(chan_idx);
            b_sample = b_all(chan_idx);
            c_sample = c_all(chan_idx);
            std_sample = std_all(chan_idx);
            I_sample = total_charge(c)/num_elecs(e);

            RT_all = a_sample.*exp(-b_sample.*I_sample)+c_sample +randn(size(std_sample)).*std_sample;

            % store min RT for each combo
            RT_out(e,c,:) = min(RT_all,[],2);
        end
    end

% plot RT_out data
    mean_RT = mean(RT_out,3);
    std_err_RT = std(RT_out,[],3)/sqrt(num_runs_per_condition);
    colors = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]/255;
    
    f=figure();
    hold on;
    for c=1:numel(total_charge)
        plot(num_elecs,squeeze(mean_RT(:,c,:)),'.','markersize',24,'color',colors(c,:))
    end
    for c=1:numel(total_charge)
        plot(repmat(num_elecs',1,2)',(squeeze(mean_RT(:,c,:))+squeeze(std_err_RT(:,c,:)).*[-1,1])','linewidth',1.5,'color',colors(c,:))
    end


