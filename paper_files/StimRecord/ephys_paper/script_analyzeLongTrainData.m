%% get experimental long train data

    exp_input_data.home_computer = 1;
    
    [exp_amp_freq_data,exp_intermittent_data] = getExperimentLongTrainData(exp_input_data);
    
    exp_amp_freq_data.stim_chan = adjustArrayDataSpikeTimes(exp_amp_freq_data.stim_chan, 0.453/1000); % stim pulse length
    exp_amp_freq_data.nonstim_chan = adjustArrayDataSpikeTimes(exp_amp_freq_data.nonstim_chan, 0.453/1000); % stim pulse length
    
    exp_intermittent_data.high_freq.stim_chan = adjustArrayDataSpikeTimes(exp_intermittent_data.high_freq.stim_chan, 0.453/1000); % stim pulse length
    exp_intermittent_data.high_freq.nonstim_chan = adjustArrayDataSpikeTimes(exp_intermittent_data.high_freq.nonstim_chan, 0.453/1000); % stim pulse length
    
    exp_intermittent_data.low_freq.stim_chan = adjustArrayDataSpikeTimes(exp_intermittent_data.low_freq.stim_chan, 0.453/1000); % stim pulse length
    exp_intermittent_data.low_freq.nonstim_chan = adjustArrayDataSpikeTimes(exp_intermittent_data.low_freq.nonstim_chan, 0.453/1000); % stim pulse length


%% plot amp-freq PSTH for each condition
    for u = 1:numel(exp_amp_freq_data.stim_chan)
        input_data.amp_freq = 1;

        input_data.window = [-1500,8000];
        input_data.unit_idx = u;
        input_data.num_cols = 3;
        input_data.account_for_artifact = 1;
        plotPSTHArrayData(exp_amp_freq_data.stim_chan,input_data);
%         f=figure(1);
%         saveFiguresLIB(f,fpath,f.Name);
%         close all
        
    end
  
%% plot intermittent_180Hz PSTH for each condition
    for u = 14%1:numel(exp_intermittent_data.high_freq.stim_chan)
        input_data.amp_freq = 0;

        input_data.window = [-1500,8000];
        input_data.unit_idx = u;
        input_data.num_cols = 3;
        input_data.account_for_artifact = 1;
        exp_intermittent_data.high_freq.stim_chan{u}.num_stims = exp_intermittent_data.high_freq.stim_chan{u}.numStims;
        plotPSTHArrayData(exp_intermittent_data.high_freq.stim_chan,input_data);
%         f=figure(1);
%         saveFiguresLIB(f,fpath,f.Name);
%         close all
    end
    
%% Get all relevant data (decay time constant and response amp) for amp-freq and intermittent data. This can take awhile....
    analyze_amp_freq_data = 1;
    analyze_intermittent_data = 0;

    decay_rate_input_data = [];
    
    decay_rate_input_data.bin_size = 50; % ms
    decay_rate_input_data.min_rate = 0.5; % Hz
    
    decay_rate_input_data.response_amp_time = 3500; % ms, ignored if num_pulses > 0
    decay_rate_input_data.response_amp_num_pulses = -1; % set as a positive number to override response_amp_time
    decay_rate_input_data.response_amp_pulse_window = [0.5,5]; % if using num_pulses, this determines when after each pulse to count spikes
    
    field_names = {'stim_chan','nonstim_chan'};
    for i_field = 1:numel(field_names)
        if(analyze_amp_freq_data)
            decay_rate_input_data.is_intermittent = 0;
            exp_amp_freq_data = getDecayRateWrapper(exp_amp_freq_data,field_names{i_field},decay_rate_input_data);
        end
        if(analyze_intermittent_data)
            decay_rate_input_data.is_intermittent = 1;
            exp_intermittent_data.high_freq = getDecayRateWrapper(exp_intermittent_data.high_freq,field_names{i_field},decay_rate_input_data);

            decay_rate_input_data.is_intermittent = 1;
            exp_intermittent_data.low_freq = getDecayRateWrapper(exp_intermittent_data.low_freq,field_names{i_field},decay_rate_input_data);
        end
    end
    
%% plot decay time constant for amp freq data for stim and non stim channels
    amp_freq_data = 1;
    marker_size = 8;
    offset = [0.95,1,1.05];
        
    colors = inferno(4);
    
    f=figure('Position',[2426 347 445 330]); hold on
    f.Name = 'stim_channel_ampfreq_decay_rate';
    
    boxplot_params.linewidth = 1.75;
    boxplot_params.box_width = 3.5;
    boxplot_params.whisker_width = boxplot_params.box_width*0.2;
    boxplot_params.outlier_marker = 'x';
    boxplot_params.outlier_marker_size = 5;
    boxplot_params.use_log_x_scale = 0;

    x_data = [131,131,131,104,104,104,80,80,80,51,51,51]+[-3.5,0,3.5,-3.5,0,3.5,-3.5,0,3.5,-3.5,0,3.5];
    color_idx = [0,1,2,0,1,2,0,1,2,0,1,2]+1;
    for condition = 1:12
        boxplot_params.outlier_color = colors(color_idx(condition),:);
        boxplot_params.median_color = colors(color_idx(condition),:);
        boxplot_params.box_color = colors(color_idx(condition),:);
        boxplot_params.whisker_color = colors(color_idx(condition),:);
        boxplot_wrapper(x_data(condition),(exp_amp_freq_data.stim_chan_decay_rates(:,condition)),boxplot_params);
    end
    xlabel('Frequency (Hz)');
    ylabel('Decay rate (1/s)')

    xlim([25,155])
    ylim([-0.5,10])
    ax = gca;
    ax.XTick = [51,80,104,131];
    formatForLee(gcf);
    set(gca,'fontsize',14);
    ax.XMinorTick = 'off';
        
%% do statistics for amp-freq decay rate
    amp_data = [20,40,60,20,40,60,20,40,60,20,40,60];
    freq_data = [131,131,131,104,104,104,80,80,80,51,51,51];
    
    % make table of decay rates and inputs for fitlm
    amp_list = []; freq_list = []; monkey_list = []; is_stim_chan = []; decay_list = [];
    unit_id = []; chan_stim = []; chan_rec = [];
    for i_unit = 1:size(exp_amp_freq_data.stim_chan_decay_rates,1)
        decay_list = [decay_list; exp_amp_freq_data.stim_chan_decay_rates(i_unit,:)'];
        amp_list = [amp_list;amp_data'];
        freq_list = [freq_list;freq_data'];
        monkey_list = [monkey_list; exp_amp_freq_data.stim_chan_monkey(i_unit)*ones(12,1)];
        is_stim_chan = [is_stim_chan; ones(12,1)];
        unit_id = [unit_id; i_unit*ones(12,1)];
        chan_stim = [chan_stim; exp_amp_freq_data.stim_chan_chan_stim(i_unit)*ones(12,1)];
        chan_rec = [chan_rec; exp_amp_freq_data.stim_chan_chan_rec(i_unit)*ones(12,1)];
    end
    for i_unit = 1:size(exp_amp_freq_data.nonstim_chan_decay_rates,1)
        decay_list = [decay_list; exp_amp_freq_data.nonstim_chan_decay_rates(i_unit,:)'];
        amp_list = [amp_list;amp_data'];
        freq_list = [freq_list;freq_data'];
        monkey_list = [monkey_list; exp_amp_freq_data.nonstim_chan_monkey(i_unit)*ones(12,1)];
        is_stim_chan = [is_stim_chan; zeros(12,1)];
        unit_id = [unit_id; (i_unit+100)*ones(12,1)]; % + 100 to avoid the stim chan id's
        chan_stim = [chan_stim; exp_amp_freq_data.nonstim_chan_chan_stim(i_unit)*ones(12,1)];
        chan_rec = [chan_rec; exp_amp_freq_data.nonstim_chan_chan_rec(i_unit)*ones(12,1)];
    end
    
    keep_mask = ones(length(decay_list),1);
    keep_mask(decay_list > 10) = 0;
    keep_mask(isnan(decay_list)) = 0;
    decay_list = decay_list(keep_mask==1);
    amp_list = amp_list(keep_mask==1);
    freq_list = freq_list(keep_mask==1);
    monkey_list = monkey_list(keep_mask==1);
    is_stim_chan = is_stim_chan(keep_mask==1);
    unit_id = unit_id(keep_mask==1);
    chan_stim = chan_stim(keep_mask==1);
    chan_rec = chan_rec(keep_mask==1);
    
    
    % compare decay rates for stim and non stim channel (with anova or with
    % a different test?)
    amp_freq_decay_tbl = table((decay_list),amp_list,freq_list,categorical(monkey_list),categorical(is_stim_chan),categorical(unit_id),...
        categorical(chan_stim),categorical(chan_rec),...
        'VariableNames',{'decay','amp','freq','monkey','is_stim_chan','unit_id','chan_stim','chan_rec'});
    
    mdl_spec = 'decay~amp*freq*is_stim_chan +monkey+chan_stim+chan_rec';
%     mdl_spec = 'decay~amp*freq + monkey + chan_stim';
    amp_freq_decay_mdl = fitlm(amp_freq_decay_tbl,mdl_spec);
    
    
    str_find = {'amp','freq','is_stim_chan','monkey_1'};
    keep_mask = zeros(size(amp_freq_decay_mdl.CoefficientNames));
    for i_str = 1:numel(str_find)
        idx_keep = find(~cellfun(@isempty,strfind(amp_freq_decay_mdl.CoefficientNames,str_find{i_str})));
        keep_mask(idx_keep)=1;
    end
    
    disp(amp_freq_decay_mdl.Formula)
    amp_freq_decay_mdl.Coefficients(keep_mask==1,:)
    disp(amp_freq_decay_mdl.Rsquared)
    
%% do rank sum tests for each stimulation condition
    unique_amps = unique(amp_list);
    unique_freqs = unique(freq_list);
    alpha = 0.05/(numel(unique_amps)*numel(unique_freqs));
    
    p_vals = zeros(numel(unique_amps),numel(unique_freqs));
    
    for i_amp = 1:numel(unique_amps)
        for i_freq = 1:numel(unique_freqs)
            stim_chan_mask = amp_list==unique_amps(i_amp) & freq_list == unique_freqs(i_freq) & is_stim_chan == 1;
            nonstim_chan_mask = amp_list==unique_amps(i_amp) & freq_list == unique_freqs(i_freq) & is_stim_chan == 0;
            
            p_vals(i_amp,i_freq) = ranksum(decay_list(stim_chan_mask), decay_list(nonstim_chan_mask),'tail','right');
        end
    end
    
%% plot decay rate for intermittent data
    boxplot_params.linewidth = 1.75;
    boxplot_params.box_width = 1.05;
    boxplot_params.whisker_width = boxplot_params.box_width-0.02;
    boxplot_params.outlier_marker = '.';
    boxplot_params.outlier_marker_size = 12;
    boxplot_params.use_log_x_scale = 1;

    x_data = [50,50,50,100,100,100,200,200,200,4000,4000,4000].*[1.1,1,0.9,1.1,1,0.9,1.1,1,0.9,1.1,1,0.9];
    color_idx = [2,1,0,2,1,0,2,1,0,2,1,0];
    f=figure('Position',[573 477 890 350]); hold on
    f.Name = 'stim_channel_intermittent_decay_rate_131Hz';
    ax1=subplot(1,2,1);
    for condition = 1:9
        boxplot_params.outlier_color = getColorFromList(1,color_idx(condition));
        boxplot_params.median_color = getColorFromList(1,color_idx(condition));
        boxplot_params.box_color = getColorFromList(1,color_idx(condition));
        boxplot_params.whisker_color = getColorFromList(1,color_idx(condition));
        boxplot_wrapper(x_data(condition),exp_intermittent_data.low_freq.stim_chan_decay_rates(:,condition),boxplot_params);
    end
    xlabel('Pulse active time (ms)');
    ylabel('Decay rate (1/s)')
    set(gca,'XScale','log')
    xlim([35,270])
    ax = gca;
    ax.XTick = [50,100,200];
    formatForLee(gcf);
    set(gca,'fontsize',14);
    ax.XMinorTick = 'off';

    ax2=subplot(1,2,2); hold on
    for condition = 10:12
        boxplot_params.outlier_color = getColorFromList(1,color_idx(condition));
        boxplot_params.median_color = getColorFromList(1,color_idx(condition));
        boxplot_params.box_color = getColorFromList(1,color_idx(condition));
        boxplot_params.whisker_color = getColorFromList(1,color_idx(condition));
        boxplot_wrapper(x_data(condition),exp_intermittent_data.low_freq.stim_chan_decay_rates(:,condition),boxplot_params);
    end
    xlabel('Pulse active time (ms)');
    ylabel('Decay rate (1/s)')
    set(gca,'XScale','log')
    xlim([3000,2.3143E4])
    ax = gca;
    ax.XTick = [4000];
    formatForLee(gcf);
    set(gca,'fontsize',14);
    ax.XMinorTick = 'off';
    linkaxes([ax1,ax2],'y');
    ax1.YLim=[-0.5,4];
    
%% do statistics for intermittent decay rate
    duty_data = [67,50,33,67,50,33,67,50,33,67,50,33];
    dur_data = [50,50,50,100,100,100,200,200,200,-1,-1,-1];
    
    % make table of decay rates and inputs for fitlm
    duty_list = []; dur_list = []; monkey_list = []; decay_list = []; freq_list = [];
    is_cont_list = []; chan_stim = [];
    for i_unit = 1:size(exp_intermittent_data.high_freq.stim_chan_decay_rates,1)
        decay_list = [decay_list; exp_intermittent_data.high_freq.stim_chan_decay_rates(i_unit,:)'];
        duty_list = [duty_list;duty_data'];
        dur_list = [dur_list;dur_data'];
        monkey_list = [monkey_list; exp_intermittent_data.high_freq.stim_chan_monkey(i_unit)*ones(12,1)];
        freq_list = [freq_list; 179*ones(12,1)];
        is_cont_list = [is_cont_list; zeros(9,1); ones(3,1)];
        chan_stim = [chan_stim; exp_intermittent_data.high_freq.stim_chan_chan_stim(i_unit)*ones(12,1)];
    end
    for i_unit = 1:size(exp_intermittent_data.low_freq.stim_chan_decay_rates,1)
        decay_list = [decay_list; exp_intermittent_data.low_freq.stim_chan_decay_rates(i_unit,:)'];
        duty_list = [duty_list;duty_data'];
        dur_list = [dur_list;dur_data'];
        monkey_list = [monkey_list; exp_intermittent_data.low_freq.stim_chan_monkey(i_unit)*ones(12,1)];
        freq_list = [freq_list; 131*ones(12,1)];
        is_cont_list = [is_cont_list; zeros(9,1); ones(3,1)];
        chan_stim = [chan_stim; exp_intermittent_data.low_freq.stim_chan_chan_stim(i_unit)*ones(12,1)];
    end
    
    keep_mask = ~isnan(decay_list);
    decay_list = decay_list(keep_mask==1);
    duty_list = duty_list(keep_mask==1);
    dur_list = dur_list(keep_mask);
    monkey_list = monkey_list(keep_mask);
    freq_list = freq_list(keep_mask);
    is_cont_list = is_cont_list(keep_mask);
    chan_stim = chan_stim(keep_mask);
    
% make two tables, one for studying how duration and duty cycle affect the data
    keep_mask = is_cont_list == 0;
    duty_decay_tbl = table(decay_list(keep_mask),duty_list(keep_mask),dur_list(keep_mask),...
        categorical(monkey_list(keep_mask)),categorical(freq_list(keep_mask)),categorical(chan_stim(keep_mask)),...
        'VariableNames',{'decay','duty','dur','monkey','freq','chan_stim'});
    
    mdl_spec = 'decay~duty + dur + monkey + freq + chan_stim';
    duty_decay_mdl = fitlm(duty_decay_tbl,mdl_spec)
    
%     % second is for comparing decay rate during intermittent and continuous
%     % stimulation
    int_decay_tbl = table(decay_list,duty_list,dur_list,...
        categorical(monkey_list),categorical(freq_list),categorical(is_cont_list),categorical(chan_stim),...
        'VariableNames',{'decay','duty','dur','monkey','freq','is_cont','chan_stim'});
    
    mdl_spec = 'decay~duty  + is_cont + monkey + freq + chan_stim';
    int_decay_mdl = fitlm(int_decay_tbl,mdl_spec)
   
    

%% run a wilcoxon rank-sum test for each stimulation frequency and decay rate (combine durations)
    unique_duty = unique(duty_list);
    unique_freq = unique(freq_list);
    
    p_vals = ones(numel(unique_duty),numel(unique_freq));
    
    for i_duty = 1:numel(unique_duty)
        for i_freq = 1:numel(unique_freq)
            int_mask = duty_list==unique_duty(i_duty) & freq_list==unique_freq(i_freq) & is_cont_list == 0;
            cont_mask = duty_list==unique_duty(i_duty) & freq_list==unique_freq(i_freq) & is_cont_list == 1;
            
            p_vals(i_duty,i_freq) = ranksum(decay_list(int_mask), decay_list(cont_mask),'tail','both');
        end
    end
    
    


    
    
%% plot response_amp vs. distance and fit 
    f=figure();
    f.Name = 'AmpFreq_response_amplitude_distance';    
    
    plot_all_conditions = 0;
    
    condition_map = [10,11,12,7,8,9,4,5,6,1,2,3]; % reorder conditions since the data order is weird
    if(plot_all_conditions)
        idx_plot = 1:12;
        subplot_map = condition_map;
        subplot_size = [4,3];
    else
%         idx_plot = [1,3,10,12]; 
        idx_plot = [4,6,10,12];
        subplot_map = [3,4,1,2];
        subplot_size = [2,2];
    end  
    
    amp_freq_fit = {};
    a_params = zeros(4,3);
    a_bounds = zeros(4,3,2);
    b_params = zeros(4,3);
    b_bounds = zeros(4,3,2);
    
    x_pred = [0:10:4500];
    for i = 1:12 % get fit and plot dot plot of responses for idx_plot
        % fit data
        amp_freq_fit{i} = fit(exp_amp_freq_data.nonstim_chan_distance_from_stim,exp_amp_freq_data.nonstim_chan_response_amp(:,i),...
            'a*exp(-x/b)','StartPoint',[50,4000],'upper',[500,1E5]);
        a_params(ceil(condition_map(i)/3),mod(condition_map(i)-1,3)+1) = amp_freq_fit{i}.a;
        b_params(ceil(condition_map(i)/3),mod(condition_map(i)-1,3)+1) = amp_freq_fit{i}.b;
        temp = confint(amp_freq_fit{i});
        a_bounds(ceil(condition_map(i)/3),mod(condition_map(i)-1,3)+1,:) = temp(:,1);
        b_bounds(ceil(condition_map(i)/3),mod(condition_map(i)-1,3)+1,:) = temp(:,2);
        
        if(any(i==idx_plot))
            % plot
            data_idx = find(i==idx_plot);
            ax(data_idx) = subplot(subplot_size(1),subplot_size(2),subplot_map(data_idx));
            plot(exp_amp_freq_data.nonstim_chan_distance_from_stim,exp_amp_freq_data.nonstim_chan_response_amp(:,i),'k.')
            hold on
            plot(x_pred,feval(amp_freq_fit{i},[0:10:4500]),'r--','linewidth',2)
            
            formatForLee(gcf);
            set(gca,'fontsize',14)
            if(i==idx_plot(1))
                xlabel('Distance (\mum)')
                ylabel('FR above baseline (Hz)');
            end
        end
    end

    linkaxes(ax,'xy');
    ylim([-30,100])
    xlim([0,4700])
    

%% new amp vs distance metrics
% (a) compare change in response amp in activated neurons as params increase
    f=figure('Position',[680 558 356 420]);
    f.Name = 'Han_Duncan_nonstimchan_deltaFR_stimparams';
    
    condition_map = [10,11,12,7,8,9,4,5,6,1,2,3]; % reorder conditions since the data order is weird
    amp_data = [20,40,60,20,40,60,20,40,60,20,40,60];
    freq_data = [131,131,131,104,104,104,80,80,80,51,51,51];
    min_cond_idx = 10;
    unique_amps = unique(amp_data);
    color_list = inferno(4);

    boxplot_params = [];
    boxplot_params.box_width = 2.35*2;
    boxplot_params.linewidth = 2;
    
    activated_mask = exp_amp_freq_data.nonstim_chan_is_responsive(:,min_cond_idx) == 1;% | exp_amp_freq_data.nonstim_chan_is_responsive(:,min_cond_idx) == 0;
    resp_data = exp_amp_freq_data.nonstim_chan_response_amp(activated_mask==1,:);
    std_resp_data = exp_amp_freq_data.nonstim_chan_std_response_amp(activated_mask==1,:);
    
    offset = [-2.6,0,2.6]*2;
    
    % data matrices for stat table
    resp_amp_all = [];
    chan_rec_all = [];
    chan_stim_all = [];
    monkey_all = [];
    amp_all = []; freq_all = [];
    for condition = 1:numel(amp_data)
            color_idx = find(unique_amps==amp_data(condition));
            pair_change = std_resp_data(:,condition).^2;%/freq_data(condition); % - resp_data(:,min_cond_idx);

            boxplot_params.outlier_color = color_list(color_idx,:);
            boxplot_params.median_color = color_list(color_idx,:);
            boxplot_params.box_color = color_list(color_idx,:);
            boxplot_params.whisker_color = color_list(color_idx,:);

            boxplot_wrapper(freq_data(condition)+offset(color_idx),pair_change,boxplot_params);
    %         errorbar(freq_data(condition)+offset(find(unique_amps==amp_data(condition))),mean(pair_change),std(pair_change),...
    %             'color',color_list(find(unique_amps==amp_data(condition)),:),'marker','.','markersize',24);
            hold on;

            % put data in stat matrices
            resp_amp_all = [resp_amp_all; std_resp_data(:,condition)];%/freq_data(condition)];
            monkey_all = [monkey_all; exp_amp_freq_data.nonstim_chan_monkey(activated_mask==1)];
            chan_stim_all = [chan_stim_all; exp_amp_freq_data.nonstim_chan_chan_stim(activated_mask==1) + 100*(exp_amp_freq_data.nonstim_chan_monkey(activated_mask==1)==1)];
            chan_rec_all = [chan_rec_all; exp_amp_freq_data.nonstim_chan_chan_rec(activated_mask==1) + 100*(exp_amp_freq_data.nonstim_chan_monkey(activated_mask==1)==1)];
            amp_all = [amp_all; amp_data(condition)*ones(size(resp_data,1),1)];
            freq_all = [freq_all; freq_data(condition)*ones(size(resp_data,1),1)];
    end
        
    xlabel('Frequency (Hz)');
    ylabel('Firing rate variance across trains (Hz^2)');
    
    resp_amp_all(resp_amp_all < 0) = 0;
    
    formatForLee(gcf);
    set(gca,'fontsize',14);
    xlim([40,142]);
    set(gca,'XTick',unique(freq_data),'XMinorTick','off');
    
    data_table = table((resp_amp_all),categorical(chan_rec_all),categorical(chan_stim_all),categorical(monkey_all),amp_all,freq_all,'VariableNames',...
        {'resp','chan_rec','chan_stim','monkey','amp','freq'});
    mdlspec = 'resp~amp+freq + chan_rec*chan_stim + monkey';
    
    resp_amp_mdl = fitlm(data_table,mdlspec);
    str_find = {'amp','freq','dist','monkey'};
    keep_mask = zeros(size(resp_amp_mdl.CoefficientNames));
    for i_str = 1:numel(str_find)
        idx_keep = find(~cellfun(@isempty,strfind(resp_amp_mdl.CoefficientNames,str_find{i_str})));
        keep_mask(idx_keep)=1;
    end
    
    disp(resp_amp_mdl.Formula)
    resp_amp_mdl.Coefficients(keep_mask==1,:)
    disp(resp_amp_mdl.Rsquared)
    
%% (a2) compare effect of amplitude for each neuron individually, pooling across frequencies.
    n = numel(exp_amp_freq_data.nonstim_chan_response_to_each_pulse);
    p_vals = nan(n,1);
    amp_data = [20,40,60,20,40,60,20,40,60,20,40,60];
    freq_data = [131,131,131,104,104,104,80,80,80,51,51,51];
    low_mask = amp_data==20;
    high_mask = amp_data==60;
    low_idx = 10;
    alpha = 0.001;
    n_pulses = 20;
    is_responsive=exp_amp_freq_data.nonstim_chan_is_responsive(:,low_idx);
    
    for i_unit = 1:n
        if(is_responsive(i_unit)==1 || 1==1)
            x=exp_amp_freq_data.nonstim_chan_response_to_each_pulse{i_unit}(high_mask);
            x=cellfun(@(y) y(1:n_pulses), x, 'un', 0);
            high_resp = vertcat(x{:}); 
            x=exp_amp_freq_data.nonstim_chan_response_to_each_pulse{i_unit}(low_mask);
            x=cellfun(@(y) y(1:n_pulses), x, 'un', 0);
            low_resp = vertcat(x{:});
% %             high_resp = exp_amp_freq_data.nonstim_chan_response_amp(i_unit,high_mask);
% %             low_resp = exp_amp_freq_data.nonstim_chan_response_amp(i_unit,low_mask);
% 
            p_vals(i_unit) = ranksum(high_resp,low_resp,'tail','right');
        end
    end
    
    sum(p_vals<alpha)
    
%% (b) compare number of activated neurons for each condition
    figure();
    offset_width = 5;
    offset = repmat([-1,0,1],1,4)*offset_width;
    num_responsive = sum(exp_amp_freq_data.nonstim_chan_is_responsive_nonstim);
    total_neurons = size(exp_amp_freq_data.nonstim_chan_is_responsive_nonstim,1);
    for condition = 1:size(exp_amp_freq_data.nonstim_chan_response_amp,2)
%         plot(freq_data(condition),num_responsive(condition),'.',...
%             'markersize',20,'color',color_list(find(amp_data(condition)==unique_amps),:));
        bar(freq_data(condition)+offset(condition),num_responsive(condition)/total_neurons,offset_width,...
            'FaceColor',color_list(find(amp_data(condition)==unique_amps),:))
        hold on;
        
    end
    
    ax=gca;
    ax.XTick = unique(freq_data);
    formatForLee(gcf);
    ax.XMinorTick = 'off';
    xlabel('Frequency (Hz)');
    ylabel('Proportion of responsive neurons');
    set(ax,'fontsize',14);
    xlim([min(freq_data)-offset_width*3,max(freq_data)+offset_width*3]);
    
    data_table = table(num_responsive',amp_data',freq_data','variableNames',{'perc','amp','freq'});
    mdlspec = 'perc~amp+freq';
    
    perc_mdl = fitglm(data_table,mdlspec,'Distribution','Binomial','BinomialSize',total_neurons)
    
%% (c) proportion of neurons activated at each distance for each param
    f=figure('Position',[680 558 400 420]);
    f.Name = 'Han_Duncan_nonstimchan_propact_dist';
    
    % data matrices for stat table
    bin_size = 400;
    bin_edges = [200:bin_size:4000];
    dist_all = [];
    chan_rec_all = [];
    chan_stim_all = [];
    monkey_all = [];
    amp_all = []; freq_all = []; prop_all = []; num_tot_all = [];
    subplot_idx = [1,2,3,1,2,3,1,2,3,1,2,3];
    color_list = inferno(4);
    unique_freq = unique(freq_data);
    for condition = 1:numel(amp_data)
        activated_mask = exp_amp_freq_data.nonstim_chan_is_responsive(:,condition) == 1;
        num_act = histcounts(exp_amp_freq_data.nonstim_chan_distance_from_stim(activated_mask),bin_edges);
        num_tot = histcounts(exp_amp_freq_data.nonstim_chan_distance_from_stim,bin_edges);
        
        prop_activated = num_act./num_tot;
        
        if(freq_data(condition) == 51 || freq_data(condition) == 131 || freq_data(condition) == 104)
            if(freq_data(condition)==51)
                ls = '-';
            elseif(freq_data(condition)==104)
                ls = '--';
            else
                ls =':';
            end
            plot(bin_edges(1:end-1)+mode(diff(bin_edges)), smooth(prop_activated,'lowess'),...
                'color',color_list(find(unique_amps==amp_data(condition)),:),'linestyle',ls,'linewidth',2)
            hold on;
        end
        
        % do stats on non-smoothed data
        num_act = histcounts(exp_amp_freq_data.nonstim_chan_distance_from_stim(activated_mask),bin_edges);
        num_tot = histcounts(exp_amp_freq_data.nonstim_chan_distance_from_stim,bin_edges);
        
        prop_activated = num_act./num_tot;
        
        % put data in stat matrices
        dist_all = [dist_all; bin_edges(1:end-1)'+mode(diff(bin_edges))];
        prop_all = [prop_all; prop_activated'];
        amp_all = [amp_all; amp_data(condition)*ones(numel(prop_activated),1)];
        freq_all = [freq_all; freq_data(condition)*ones(numel(prop_activated),1)];
        num_tot_all = [num_tot_all; num_tot'];
    end
    
    xlabel('Distance (\mum)');
    ylabel('Proportion of activated neurons');
    
    formatForLee(gcf);
    set(gca,'fontsize',14);
    
    data_table = table(dist_all,prop_all,freq_all,amp_all,...
        'VariableNames',{'dist','prop','freq','amp'});
    mdlspec = 'prop~dist + amp + freq';
    
    data_table.prop = data_table.prop.*num_tot_all;
    dist_mdl = fitglm(data_table,mdlspec,'Distribution','Binomial','BinomialSize',num_tot_all)
%     dist_mdl = stepwiseglm(data_table,mdlspec,'Distribution','Binomial','BinomialSize',num_tot_all)
%     dist_mdl = fitlm(data_table,mdlspec)
%     str_find = {'amp','freq','dist','monkey_1','(Intercept)'};
%     keep_mask = zeros(size(dist_mdl.CoefficientNames));
%     for i_str = 1:numel(str_find)
%         idx_keep = find(~cellfun(@isempty,strfind(dist_mdl.CoefficientNames,str_find{i_str})));
%         keep_mask(idx_keep)=1;
%     end
%     
%     disp(dist_mdl.Formula)
%     dist_mdl.Coefficients(keep_mask==1,:)
%     disp(dist_mdl.Rsquared)
    
    
%% (d) response amp with distance for activated neurons
    figure();
    % data matrices for stat table
    bin_size = 400;
    bin_edges = [200:bin_size:4000];
    dist_all = [];
    chan_rec_all = [];
    chan_stim_all = [];
    monkey_all = [];
    amp_all = []; freq_all = []; prop_all = [];
    subplot_idx = [1,2,3,1,2,3,1,2,3,1,2,3];
    color_list = inferno(4);
    unique_freq = unique(freq_data);
    for condition = 1:numel(amp_data)
        activated_mask = exp_amp_freq_data.nonstim_chan_is_responsive(:,condition) == 1;
        resp_amp = exp_amp_freq_data.nonstim_chan_is_responsive(activated_mask,:);
        dist = exp_amp_freq_data.nonstim_chan_distance_from_stim(activated_mask);
        
        if(freq_data(condition) == 51 || freq_data(condition) == 104)
            if(freq_data(condition)==51)
                ls = '-';
            else
                ls ='--';
            end
            plot(bin_edges(1:end-1)+mode(diff(bin_edges)), prop_activated,...
                'color',color_list(find(unique_amps==amp_data(condition)),:),'linestyle',ls,'linewidth',2)
            hold on;
        end
        % put data in stat matrices
        dist_all = [dist_all; bin_edges(1:end-1)'+mode(diff(bin_edges))];
        prop_all = [prop_all; prop_activated'];
        amp_all = [amp_all; amp_data(condition)*ones(numel(prop_activated),1)];
        freq_all = [freq_all; freq_data(condition)*ones(numel(prop_activated),1)];
    end
    
    xlabel('Distance (\mum)');
    ylabel('Proportion of activated neurons');
    
    formatForLee(gcf);
    set(gca,'fontsize',14);
    
    data_table = table(dist_all,prop_all,freq_all,amp_all,...
        'VariableNames',{'dist','prop','freq','amp'});
    mdlspec = 'prop~dist*amp+dist*freq';
    
    dist_mdl = fitlm(data_table,mdlspec);
    
    str_find = {'amp','freq','dist','monkey_1','(Intercept)'};
    keep_mask = zeros(size(dist_mdl.CoefficientNames));
    for i_str = 1:numel(str_find)
        idx_keep = find(~cellfun(@isempty,strfind(dist_mdl.CoefficientNames,str_find{i_str})));
        keep_mask(idx_keep)=1;
    end
    
    disp(dist_mdl.Formula)
    dist_mdl.Coefficients(keep_mask==1,:)
    disp(dist_mdl.Rsquared)
    
    


%% compare response amp vs distance metrics. 
% Fit with exponential decay.
% Longest distance where x percentile (75th?) response amp is above some value?

%     f=figure();
%     f.Name = 'AmpFreq_response_amplitude_distance';    
%     
    plot_all_conditions = 0;
    
    condition_map = [10,11,12,7,8,9,4,5,6,1,2,3]; % reorder conditions since the data order is weird
    if(plot_all_conditions)
        idx_plot = 1:12;
        subplot_map = condition_map;
        subplot_size = [4,3];
    else
        idx_plot = [4,6,10,12]; % 20uA 104Hz, 60uA 104Hz, 20uA 51Hz, 60uA 51Hz.
        subplot_map = [3,4,1,2];
        subplot_size = [2,2];
    end  
    
    amp_freq_fit = {};
    a_params = zeros(4,3);
    a_bounds = zeros(4,3,2);
    b_params = zeros(4,3);
    b_bounds = zeros(4,3,2);
    last_bin = zeros(4,3);
    
    prctile_val = 75;
    response_amp_threshold = 30;
    
    x_pred = [0:10:4500];
    x_step = 100;
    x_width = 200; % 200 both ways
    x_max = 4500;
    for i = 1:12 % get fit and plot dot plot of responses for idx_plot
        % fit data
        amp_freq_fit{i} = fit(exp_amp_freq_data.nonstim_chan_distance_from_stim,exp_amp_freq_data.nonstim_chan_response_amp(:,i),...
            'a*exp(-x/b)','StartPoint',[50,4000],'upper',[500,1E5]);
        a_params(ceil(condition_map(i)/3),mod(condition_map(i)-1,3)+1) = amp_freq_fit{i}.a;
        b_params(ceil(condition_map(i)/3),mod(condition_map(i)-1,3)+1) = amp_freq_fit{i}.b;
        temp = confint(amp_freq_fit{i});
        a_bounds(ceil(condition_map(i)/3),mod(condition_map(i)-1,3)+1,:) = temp(:,1);
        b_bounds(ceil(condition_map(i)/3),mod(condition_map(i)-1,3)+1,:) = temp(:,2);
        
        % get largest distance where x percentile response amp is above
        % some value
        % take steps
        x_curr = x_width;
        prctile_data = []; x_data = [];
        while(x_curr < x_max)
            resp_mask = exp_amp_freq_data.nonstim_chan_distance_from_stim > x_curr-x_width & exp_amp_freq_data.nonstim_chan_distance_from_stim < x_curr+x_width;
            response_data = exp_amp_freq_data.nonstim_chan_response_amp(resp_mask,i);
            prctile_data(end+1) = prctile(response_data,prctile_val);
            x_data(end+1) = x_curr;
            
            x_curr = x_curr+x_step;
            
        end
        if(~isempty(find(prctile_data > response_amp_threshold,1,'last')))
            last_bin(ceil(condition_map(i)/3),mod(condition_map(i)-1,3)+1) = x_data(find(prctile_data > response_amp_threshold,1,'last'));
        else
            last_bin(ceil(condition_map(i)/3),mod(condition_map(i)-1,3)+1) = 0;
        end
        
    end


    
    
    
    
    

%     % only use 1 monkey at a time?
%     for idx = unique(monkey_data)'
%         data_table = table(response_data(monkey_data==idx),distance_data(monkey_data==idx),amp_data(monkey_data==idx),...
%             freq_data(monkey_data==idx),categorical(chan_rec_data(monkey_data==idx)),categorical(chan_stim_data(monkey_data==idx)),...
%             'VariableNames',{'resp','dist','amp','freq','chan_rec','chan_stim'});
% 
%         amp_freq_mdl{end+1} = fitlm(data_table,modelspec)
% %         amp_freq_mdl{end+1} = fitglm(data_table,modelspec,'Distribution','poisson')
%     end

%     
%     for condition = 1:12
%         subplot(4,3,subplot_map(condition))
%         hold on
%         y_pred = predict(amp_freq_mdl{1},[x_pred',amps(condition)+zeros(numel(x_pred),1),freqs(condition)+zeros(numel(x_pred),1)]);
%         
%         plot(x_pred,y_pred,'k--','linewidth',2)
%     end
    
% %% resample data w/ smaller number of stimulations
% 
%     num_boot = 100;
%     num_stims_use = 8; % really should be less than what I actually used
%     
%     p_list = zeros(numel(array_data),numel(array_data{1}.binCounts),num_boot);
%     num_sig = zeros(numel(array_data),numel(array_data{1}.binCounts));
%     
%     for unit_idx = 1:numel(array_data)
%         f=figure();f.Position = [381 -83 1387 856];
%         bootstrapped_data = [];
%         bootstrapped_firing_rate = zeros(numel(array_data{unit_idx}.binCounts),num_boot,numel(array_data{unit_idx}.binCounts{1}));
% 
%         window_post_stim = [0,5];
%         pulse_data = {};
%         input_data_all{unit_idx} = input_data;
%         array_data{unit_idx}.num_stims = array_data{unit_idx}.numStims;
%         array_data{unit_idx}.trial_num = array_data{unit_idx}.stimData;
%         
%         for condition = 1:numel(array_data{unit_idx}.binCounts) % for each condition
%             if(input_data_all{unit_idx}.num_pulses(condition) == 1) % get timing of each pulse
%                 pulse_time = 0;
%             else
%                 pulse_time = [0,input_data_all{unit_idx}.IPI(condition)*(1:input_data_all{unit_idx}.num_pulses(condition)-1)];
%             end
%             pulse_data{condition}.count = zeros(num_boot,numel(pulse_time));
% 
% 
% 
%             for boot = 1:num_boot % for each bootstrap iteration
%                 % downsample data
% 
%                 stim_idx_use = datasample(1:array_data{unit_idx}.num_stims(condition),...
%                     min(num_stims_use,array_data{unit_idx}.num_stims(condition)),'Replace',false);
%                  
%                 spike_trial_times = array_data{unit_idx}.spikeTrialTimes{condition}(sum(array_data{unit_idx}.trial_num{condition} == stim_idx_use') > 0);
%                 
%                 % rebin data and store
%                 bootstrapped_firing_rate(condition,boot,:) = histcounts(spike_trial_times*1000,array_data{unit_idx}.binEdges{condition})/...
%                     (mode(diff(array_data{unit_idx}.binEdges{condition}))/1000)/num_stims_use;
%                 
%                 % get p-value for comparing early and late stim
%                 trial_counts = zeros(num_stims_use,2);
%                 for trial = 1:num_stims_use
%                     spike_trial_times = array_data{unit_idx}.spikeTrialTimes{condition}(array_data{unit_idx}.trial_num{condition} == stim_idx_use(trial)');
%                     trial_counts(trial,:) = [sum(spike_trial_times > 0 & spike_trial_times < 0.5), sum(spike_trial_times > 3 & spike_trial_times < 3.5)];
%                 end
%                 
%                 [~,p_list(unit_idx,condition,boot)] = ttest(trial_counts(:,1),trial_counts(:,2));
%                 
%                 for p = 1:numel(pulse_time)
%                     bin_idx = [find(array_data{unit_idx}.binEdges{condition} > pulse_time(p)+window_post_stim(1),1,'first'),...
%                         find(array_data{unit_idx}.binEdges{condition} > pulse_time(p)+window_post_stim(2),1,'first')];
% 
%                     pulse_data{condition}.count(boot,p) = mean(bootstrapped_firing_rate(condition,boot,bin_idx(1):bin_idx(2)));
%                 end
%             end
% 
%             subplot(3,4,condition)
%             x_data = array_data{unit_idx}.binEdges{condition}(1:end-1) + mode(diff(array_data{unit_idx}.binEdges{condition}))/2;
%             plot(x_data,squeeze(bootstrapped_firing_rate(condition,1:5:end,:)),'linewidth',1)
% 
%             subplot(3,4,condition+4)
%             x_data = array_data{unit_idx}.binEdges{condition}(1:end-1) + mode(diff(array_data{unit_idx}.binEdges{condition}))/2;
%             plot(x_data,squeeze(bootstrapped_firing_rate(condition,1,:)),'linewidth',2)
%             mean_bin_counts = squeeze(mean(bootstrapped_firing_rate(condition,:,:)));
% 
%             sorted_bootstrapped_firing_rate = sort(bootstrapped_firing_rate,2);
%             idx_use = [floor(num_boot*0.025),ceil(num_boot*0.975)];
%             firing_rate_bound = [squeeze(sorted_bootstrapped_firing_rate(condition,idx_use,:))];
%     %         
%             plot(x_data,mean_bin_counts) % plot the mean
%             hold on
%             plot(x_data, firing_rate_bound,'--')
% 
%             subplot(3,4,condition+8)
%             bin_idx = [find(array_data{unit_idx}.binEdges{condition} >= 0,1,'first'),...
%                         find(array_data{unit_idx}.binEdges{condition} >= 1000,1,'first'),...
%                         find(array_data{unit_idx}.binEdges{condition} >= 2500,1,'first'),...
%                         find(array_data{unit_idx}.binEdges{condition} >= 3500,1,'first')];
%                     
%             histogram(mean(bootstrapped_firing_rate(condition,:,bin_idx(1):bin_idx(2)),3) -...
%                 mean(bootstrapped_firing_rate(condition,:,bin_idx(3):bin_idx(4)),3));
%             
%             
%             
%         end
%     end
%        
%     num_sig = sum(p_list < 0.05,3)