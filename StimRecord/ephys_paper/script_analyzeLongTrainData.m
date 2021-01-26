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
    for u = 10:12%numel(array_data)
        input_data.amp_freq = 0;

        input_data.window = [-1500,15000];
        input_data.unit_idx = u;
        input_data.num_cols = 3;
        input_data.account_for_artifact = 1;
        plotPSTHArrayData(exp_intermittent_data.high_freq.stim_chan,input_data);
%         f=figure(1);
%         saveFiguresLIB(f,fpath,f.Name);
%         close all
    end
    
%% get decay time constant and response amp for amp-freq and intermittent data. This can take awhile....
    analyze_amp_freq_data = 0;
    analyze_intermittent_data = 1;

    decay_rate_input_data = [];
    
    decay_rate_input_data.bin_size = 50; % ms
    decay_rate_input_data.min_rate = 0.5; % Hz
    
    decay_rate_input_data.response_amp_time = 250; % ms, ignored if num_pulses > 0
    decay_rate_input_data.response_amp_num_pulses = -1; % set as a positive number to override response_amp_time
    decay_rate_input_data.response_amp_pulse_window = [1,5]; % if using num_pulses, this determines when after each pulse to count spikes
    
    field_names = {'stim_chan'};
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
    f.Name = 'nonstim_channel_ampfreq_decay_rate';
    
    boxplot_params.linewidth = 1.75;
    boxplot_params.box_width = 3.5;
    boxplot_params.whisker_width = boxplot_params.box_width*0.2;
    boxplot_params.outlier_marker = '.';
    boxplot_params.outlier_marker_size = 12;
    boxplot_params.use_log_x_scale = 0;

    x_data = [131,131,131,104,104,104,80,80,80,51,51,51]+[-3.5,0,3.5,-3.5,0,3.5,-3.5,0,3.5,-3.5,0,3.5];
    color_idx = [0,1,2,0,1,2,0,1,2,0,1,2]+1;
    for condition = 1:12
        boxplot_params.outlier_color = colors(color_idx(condition),:);
        boxplot_params.median_color = colors(color_idx(condition),:);
        boxplot_params.box_color = colors(color_idx(condition),:);
        boxplot_params.whisker_color = colors(color_idx(condition),:);
        boxplot_wrapper(x_data(condition),(exp_amp_freq_data.nonstim_chan_decay_rates(:,condition)),boxplot_params);
    end
    xlabel('Frequency (Hz)');
    ylabel('Decay rate (1/s)')

    xlim([25,155])
    ax = gca;
    ax.XTick = [51,80,104,131];
    formatForLee(gcf);
    set(gca,'fontsize',14);
    ax.XMinorTick = 'off';
        
%% plot decay rate for intermittent data
    boxplot_params.linewidth = 1.75;
    boxplot_params.box_width = 1.05;
    boxplot_params.whisker_width = boxplot_params.box_width-0.02;
    boxplot_params.outlier_marker = '.';
    boxplot_params.outlier_marker_size = 12;
    boxplot_params.use_log_x_scale = 1;

    x_data = [50,50,50,100,100,100,200,200,200,4000,4000,4000].*[1.1,1,0.9,1.1,1,0.9,1.1,1,0.9,1.1,1,0.9];
    color_idx = [2,1,0,2,1,0,2,1,0,2,1,0];
    f=figure(); hold on
    f.Name = 'stim_channel_intermittent_decay_rate_131Hz';
    f.Position = [550 473 890 420];
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
        boxplot_wrapper(x_data(condition),exp_intermittent_data.high_freq.stim_chan_decay_rates(:,condition),boxplot_params);
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
    ax1.YLim(1)=0;
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
        idx_plot = [4,6,10,12]; % 20uA 104Hz, 60uA 104Hz, 20uA 51Hz, 60uA 51Hz.
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
    ylim([-30,250])
    xlim([0,4700])
    
%% plot a (intercept) and b (space constant) across amplitudes and frequencies
    f=figure(); hold on;
    f.Name = 'Han_duncan_long_trains_spread_param_summary';
    
    offset = [-1,0,1]*2;
    
    freq_list = [51,80,104,131];
    amp_list = [20,40,60];
%     amp_colors = [0,200,0; 0,128,0; 0,100,0]/255;
    amp_colors = inferno(4);
        
    for i_plot = 1:2
        subplot(1,2,i_plot)
        switch i_plot
            case 1
                data = b_params;
                bounds = b_bounds;
            case 2
                data = a_params;
                bounds = a_bounds;
        end
        
        for i_freq = 1:numel(freq_list)
            for i_amp = 1:numel(amp_list)
                errorbar(freq_list(i_freq)+offset(i_amp),data(i_freq,i_amp),...
                    data(i_freq,i_amp)-bounds(i_freq,i_amp,1),bounds(i_freq,i_amp,2)-data(i_freq,i_amp),...
                    'marker','.','linestyle','none','color',amp_colors(i_amp,:),'markersize',18,'linewidth',1.5);
                hold on;
            end
        end
        formatForLee(gcf);
        set(gca,'fontsize',14);
        xlabel('Frequency (Hz)');
        xlim([40,140])
        
        if(i_plot == 1)
            ylabel('Space constant (1/\mum)');
            ylim([0,5000])
            l=legend('20\muA','40\muA','60\muA');
            set(l,'box','off');
        else
            ylabel('Intercept (Hz)');
            ylim([0,100])
        end
    end    
    
%% build GLM and then assess effect of amp/frequency across distances

    amps = [20,40,60,20,40,60,20,40,60,20,40,60];
    freqs = [131,131,131,104,104,104,80,80,80,51,51,51];
    
    response_data = [];
    distance_data = [];
    amp_data = [];
    freq_data = [];
    chan_rec_data = []; % add 100 if monkey is han
    monkey_data = []; % 1 = han, 0 = duncan
    chan_stim_data = []; % add 100 if monkey is han
    
    for condition = 1:numel(amps)
        distance_data = [distance_data;exp_amp_freq_data.nonstim_chan_distance_from_stim];
        response_data = [response_data;exp_amp_freq_data.nonstim_chan_response_amp(:,condition)];
        amp_data = [amp_data;zeros(numel(exp_amp_freq_data.nonstim_chan_distance_from_stim),1)+amps(condition)];
        freq_data = [freq_data;zeros(numel(exp_amp_freq_data.nonstim_chan_distance_from_stim),1)+freqs(condition)];
        monkey_data = [monkey_data; exp_amp_freq_data.nonstim_chan_monkey];
        chan_stim_data = [chan_stim_data; exp_amp_freq_data.nonstim_chan_chan_stim + 100*(exp_amp_freq_data.nonstim_chan_monkey==1)];
        chan_rec_data = [chan_rec_data; exp_amp_freq_data.nonstim_chan_chan_rec + 100*(exp_amp_freq_data.nonstim_chan_monkey==1)];
    end
    
    
    response_data(response_data < 0) = 0;
    
    modelspec = 'resp ~ dist+amp+freq+monkey+chan_rec*chan_stim';
    amp_freq_mdl = {};
    % both monkeys simultaneously
    data_table = table(response_data,distance_data,amp_data,...
            freq_data,categorical(monkey_data),categorical(chan_stim_data),categorical(chan_rec_data),...
            'VariableNames',{'resp','dist','amp','freq','monkey','chan_stim','chan_rec'});

    amp_freq_mdl{end+1} = fitlm(data_table,modelspec);
     
%% rebound excitation stats 
% duration, percent of cells
    rebound_input_data.cond_list = [1:12];
    rebound_input_data.bin_window = [-200,8000]./1000;
    rebound_input_data.pre_window = [-1000,-10]./1000;
    rebound_input_data.post_window = [10,500]./1000;
    rebound_input_data.blank_time = [-10,10]/1000; %s
    
    rebound_input_data.bin_size = 5/1000; % s
    rebound_input_data.kernel_length = 2;
    rebound_input_data.num_consec_bins = 6;
    [rebound_data] = getReboundExcitationWrapper(exp_amp_freq_data.nonstim_chan,rebound_input_data);
    
%% plot rebound excitation data -- 

    % plot duration across frequencies (with lines between same cells)
    f=figure('Position',[1403 522 395 420]);
    f.Name = 'Han_duncan_long_trains_rebound_excitation';
    suptitle('Short Trains')
    subplot(2,1,2); hold on;
    freq_list=[180,100,50,20];
    for i_cell = 1:size(rebound_data.is_rebound,1)
        plot(freq_list,rebound_data.rebound_dur(i_cell,2:end)*1000,'-k','marker','.','markersize',12);
    end
    xlim([0,200])
    ylim([0,300]);
    formatForLee(gcf);
    xlabel('Frequency (Hz)');
    ylabel('Duration (ms)');
    set(gca,'fontsize',14);
    ax=gca;
    ax.XTick = sort(freq_list);
    ax.XMinorTick = 'off';


    % also plot % of cells with rebound across frequencies
    subplot(2,1,1); hold on;
    freq_list=[180,100,50,20];
    frac_data = nan(size(freq_list));
    for i_freq = 2:size(rebound_data.is_rebound,2)
        frac_data(i_freq-1) = sum(rebound_data.is_rebound(:,i_freq),'omitnan')/sum(~isnan(rebound_data.is_rebound(:,i_freq)));
    end
    
    bar(freq_list,frac_data,'EdgeColor','k','FaceColor','k');
    xlim([0,200])
    ylim([0,1]);
    formatForLee(gcf);
    xlabel('Frequency (Hz)');
    ylabel('Fraction cells');
    set(gca,'fontsize',14);
    ax=gca;
    ax.XTick = sort(freq_list);
    ax.XMinorTick = 'off';
 
    
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
        % bin neurons by distance, get x percentile response amp in each bin
%         [~,~,bin_idx] = histcounts(exp_amp_freq_data.nonstim_chan_distance_from_stim,x_bin_edges);
%         
%         unique_bin_idx = unique(bin_idx);
%         prctile_data = nan(numel(unique_bin_idx),1);
%         for i_bin = 1:numel(unique_bin_idx)
%             response_data = exp_amp_freq_data.nonstim_chan_response_amp(bin_idx==unique_bin_idx(i_bin),i);
%             prctile_data(i_bin) = prctile(response_data,prctile_val);
%         end
%         
        if(~isempty(find(prctile_data > response_amp_threshold,1,'last')))
            last_bin(ceil(condition_map(i)/3),mod(condition_map(i)-1,3)+1) = x_data(find(prctile_data > response_amp_threshold,1,'last'));
        else
            last_bin(ceil(condition_map(i)/3),mod(condition_map(i)-1,3)+1) = 0;
        end
        
%         if(any(i==idx_plot))
%             % plot
%             data_idx = find(i==idx_plot);
%             ax(data_idx) = subplot(subplot_size(1),subplot_size(2),subplot_map(data_idx));
%             plot(exp_amp_freq_data.nonstim_chan_distance_from_stim,exp_amp_freq_data.nonstim_chan_response_amp(:,i),'k.')
%             hold on
%             plot(x_pred,feval(amp_freq_fit{i},[0:10:4500]),'r--','linewidth',2)
%             
%             formatForLee(gcf);
%             set(gca,'fontsize',14)
%             if(i==idx_plot(1))
%                 xlabel('Distance (\mum)')
%                 ylabel('FR above baseline (Hz)');
%             end
%         end
    end

%     linkaxes(ax,'xy');
%     ylim([-30,250])
%     xlim([0,4700])   
    
    
    %         amp_freq_mdl{end+1} = fitglm(data_table,modelspec,'Distribution','poisson')

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