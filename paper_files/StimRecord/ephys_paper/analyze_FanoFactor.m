%% get data    
    exp_input_data.home_computer = 1;
    
    [exp_amp_freq_data,exp_intermittent_data] = getExperimentLongTrainData(exp_input_data);
    
    exp_amp_freq_data.stim_chan = adjustArrayDataSpikeTimes(exp_amp_freq_data.stim_chan, 0.453/1000); % stim pulse length
    exp_amp_freq_data.nonstim_chan = adjustArrayDataSpikeTimes(exp_amp_freq_data.nonstim_chan, 0.453/1000); % stim pulse length
    
    
    
%% for each condition and neuron, compute a fano factor during stim and during baseline

    amp_data = [20,40,60,20,40,60,20,40,60,20,40,60];
    freq_data = [131,131,131,104,104,104,80,80,80,51,51,51];
    
    n_neurons = numel(exp_amp_freq_data.nonstim_idx);
    n_conditions = 12;
    
    mean_stim = zeros(numel(exp_amp_freq_data.nonstim_chan{1}.spikeTrialTimes) , n_neurons); % conditions across neurons
    std_stim = zeros(size(mean_stim));
    mean_base = zeros(size(mean_stim));
    std_base = zeros(size(mean_stim));
    
    for i_neuron = 1:n_neurons
        
        array_data = exp_amp_freq_data.nonstim_chan{i_neuron};
        for i_cond = 1:n_conditions
            % compute firing rate during each train of a condition
            % compute firing rate during baseline period for each condition
            n_spikes_stim = [];
            n_spikes_base = [];
            for i_trial = 1:max(array_data.stimData{i_cond})
                spike_ts = array_data.spikeTrialTimes{i_cond}(array_data.stimData{i_cond}==i_trial);
                    n_spikes_stim(end+1) = sum(spike_ts > i_w & spike_ts < i_w+1);
                    n_spikes_base(end+1) = sum(spike_ts > 10+i_w & spike_ts < 11+i_w);
            end
        
        
            % get fano factor (mean/variance) for stim and baseline
            mean_stim(i_cond,i_neuron) = mean(n_spikes_stim);
            std_stim(i_cond,i_neuron) = std(n_spikes_stim);
            mean_base(i_cond,i_neuron) = mean(n_spikes_base);
            std_base(i_cond,i_neuron) = std(n_spikes_base);
            
        end
    end
    
    
    fano_factor_stim = (std_stim.^2)./mean_stim;
    fano_factor_base = (std_base.^2)./mean_base;
    
    
%% plot distribution of fano factor for each condition
    
    offset = [0.9,1,1.1];
    colors = inferno(4);
    
    f=figure('Position',[2426 347 445 330]); hold on
%     f.Name = 'stim_channel_ampfreq_decay_rate';
    
    boxplot_params.linewidth = 1.75;
    boxplot_params.box_width = 1.5;
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
        boxplot_wrapper(x_data(condition),fano_factor_stim(condition,:),boxplot_params);
        boxplot_wrapper(x_data(condition)+1.75,fano_factor_base(condition,:), boxplot_params);
    end
        
    
    %%
    xlabel('Frequency (Hz)');
    ylabel('Decay rate (1/s)')

    xlim([25,155])
    ax = gca;
    ax.XTick = [51,80,104,131];
    formatForLee(gcf);
    set(gca,'fontsize',14);
    ax.XMinorTick = 'off';













    
    