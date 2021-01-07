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
    for u = 10:12%numel(exp_amp_freq_data.stim_chan)
        input_data.amp_freq = 1;

        input_data.window = [-1500,15000];
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
    analyze_amp_freq_data = 1;
    analyze_intermittent_data = 0;

    decay_rate_input_data = [];
    
    decay_rate_input_data.bin_size = 50; % ms
    decay_rate_input_data.min_rate = 0.5; % Hz
    
    decay_rate_input_data.response_amp_time = 250; % ms, ignored if num_pulses > 0
    decay_rate_input_data.response_amp_num_pulses = -1; % set as a positive number to override response_amp_time
    decay_rate_input_data.response_amp_pulse_window = [1,5]; % if using num_pulses, this determines when after each pulse to count spikes
    
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
        
    f=figure(); hold on
    f.Name = 'stim_channel_ampfreq_decay_rate';
    
    boxplot_params.linewidth = 1.75;
    boxplot_params.box_width = 3.5;
    boxplot_params.whisker_width = boxplot_params.box_width*0.2;
    boxplot_params.outlier_marker = '.';
    boxplot_params.outlier_marker_size = 12;
    boxplot_params.use_log_x_scale = 0;

    x_data = [131,131,131,104,104,104,80,80,80,51,51,51]+[-3.5,0,3.5,-3.5,0,3.5,-3.5,0,3.5,-3.5,0,3.5];
    color_idx = [0,1,2,0,1,2,0,1,2,0,1,2];
    f.Position = [550 473 445 420];
    for condition = 1:12
        boxplot_params.outlier_color = getColorFromList(1,color_idx(condition));
        boxplot_params.median_color = getColorFromList(1,color_idx(condition));
        boxplot_params.box_color = getColorFromList(1,color_idx(condition));
        boxplot_params.whisker_color = getColorFromList(1,color_idx(condition));
        boxplot_wrapper(x_data(condition),(exp_amp_freq_data.stim_chan_decay_rates(:,condition)),boxplot_params);
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
    f.Name = 'stim_channel_intermittent_decay_rate';
    f.Position = [550 473 890 420];
    ax1=subplot(1,2,1);
    for condition = 1:9
        boxplot_params.outlier_color = getColorFromList(1,color_idx(condition));
        boxplot_params.median_color = getColorFromList(1,color_idx(condition));
        boxplot_params.box_color = getColorFromList(1,color_idx(condition));
        boxplot_params.whisker_color = getColorFromList(1,color_idx(condition));
        boxplot_wrapper(x_data(condition),exp_intermittent_data.high_freq.stim_chan_decay_rates(:,condition),boxplot_params);
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
    
%% plot response_amp vs. distance and fit 
    f=figure();
    f.Name = 'AmpFreq_response_amplitude_distance';
    make_heatmap = 0;
    
    distance_bin_size = 500; %um, box plots and heatmap 
    x_spacing = 200;
    y_spacing = 10;
    edges = {[400-x_spacing:x_spacing:4500],[0:y_spacing:150]};
    dist_spacing = 120;
    
    idx_plot = [4,6,10,12];
    
    if(make_heatmap == 1)
%         subplot_map = [10,11,12,7,8,9,4,5,6,1,2,3];
        subplot_map = [10,11,12,3,8,4,4,5,6,1,2,2];
        x_pred = [0:10:4500];
        for i = idx_plot
            ax(i) = subplot(2,2,subplot_map(i));
%             plot(distance_from_stim,response_amp(:,i),'.')
            [N,C] = hist3([exp_amp_freq_data.nonstim_chan_distance_from_stim,exp_amp_freq_data.nonstim_chan_response_amp(:,i)],edges);

            % adjust C so it corresponds to the edge, not center
            for i_dim = 1:numel(C)
                C{i_dim} = C{i_dim} - mean(diff(C{i_dim}))/2;
            end
            pcolor(C{1},C{2},N')
            shading interp
            colormap viridis
            ax(i).CLim = [0,50];
            colorbar
            disp(sum(sum(N)))
%             [h,b]=densityplot(exp_amp_freq_data.nonstim_chan_distance_from_stim,exp_amp_freq_data.nonstim_chan_response_amp(:,i),edges)
            

            amp_freq_fit{i} = fit(exp_amp_freq_data.nonstim_chan_distance_from_stim,exp_amp_freq_data.nonstim_chan_response_amp(:,i),'a*exp(-x/b)','StartPoint',[50,4000],'upper',[500,1E5]);
            a_params(ceil(subplot_map(i)/3),mod(subplot_map(i)-1,3)+1) = amp_freq_fit{i}.a;
            b_params(ceil(subplot_map(i)/3),mod(subplot_map(i)-1,3)+1) = amp_freq_fit{i}.b;
        end

        linkaxes(ax,'xy');
%         ylim([-30,150])
        
    elseif(make_heatmap == 0) % dots for each rec:stim pair
%         subplot_map = [10,11,12,7,8,9,4,5,6,1,2,3];
        subplot_map = [10,11,12,3,8,4,4,5,6,1,2,2];
        amp_freq_fit = {};
        a_params = zeros(4,3);
        b_params = zeros(4,3);

        x_pred = [0:10:4500];
        for i = idx_plot
            ax(i) = subplot(2,2,subplot_map(i));
            plot(exp_amp_freq_data.nonstim_chan_distance_from_stim,exp_amp_freq_data.nonstim_chan_response_amp(:,i),'k.')
            amp_freq_fit{i} = fit(exp_amp_freq_data.nonstim_chan_distance_from_stim,exp_amp_freq_data.nonstim_chan_response_amp(:,i),'a*exp(-x/b)','StartPoint',[50,4000],'upper',[500,1E5]);
            hold on
            plot(x_pred,feval(amp_freq_fit{i},[0:10:4500]),'r--','linewidth',2)

            a_params(ceil(subplot_map(i)/3),mod(subplot_map(i)-1,3)+1) = amp_freq_fit{i}.a;
            b_params(ceil(subplot_map(i)/3),mod(subplot_map(i)-1,3)+1) = amp_freq_fit{i}.b;
            
            formatForLee(gcf);
            set(gca,'fontsize',14)
            if(i==idx_plot(1))
                xlabel('Distance (\mum)')
                ylabel('FR above baseline (Hz)');
            end
        end

        linkaxes(ax,'xy');
        ylim([-30,250])
        xlim([0,4700])
    end
    
    
%% build GLM

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

