%% get long train data    
    exp_input_data.home_computer = 1;
    
    [exp_amp_freq_data,exp_intermittent_data] = getExperimentLongTrainData(exp_input_data);
    
    exp_amp_freq_data.stim_chan = adjustArrayDataSpikeTimes(exp_amp_freq_data.stim_chan, 0.453/1000); % stim pulse length
    exp_amp_freq_data.nonstim_chan = adjustArrayDataSpikeTimes(exp_amp_freq_data.nonstim_chan, 0.453/1000); % stim pulse length
    
%% generate kernel for convolution
    dt = 0.0001;
    kernel_SD = 0.0005;
    % kernel half length is 3·SD out
    kernel_hl = ceil( 3 * kernel_SD / (dt) );
    % create the kernel --it will have length 2*kernel_hl+1
    kernel = normpdf( -kernel_hl*(dt) : ...
    dt : kernel_hl*(dt), ...
        0, kernel_SD );
    % compute normalization factor --this factor depends on the number of taps
    % actually used
    bin_edges = -5:dt:20;
    nm = conv(kernel,ones(1,numel(bin_edges)-1)');

%% for each stim channel (exp_amp_freq_data.nonstim_chan_idx), generate peak to peak amp 
% convolve gaussian with spikes from all neurons for each trial
% compute peak to peak amplitude after each pulse (0-5 ms), compute peak to
% peak amplitude during random baseline segments as a comparison

% this code takes 5-10 minutes to run, mostly because it's poorly
% optimized
amp_data = [20,40,60,20,40,60,20,40,60,20,40,60];
freq_data = [131,131,131,104,104,104,80,80,80,51,51,51];
ISI_edges = (0:0.1:5)/1000; % in s
unique_idx = unique(exp_amp_freq_data.nonstim_idx);
window = [1,5]/1000; % look 1-5 ms after each pulse for the maximum spike train amp

n_neurons = zeros(numel(unique_idx),1);
for i = unique_idx
    is_chan = exp_amp_freq_data.nonstim_idx == i;
    n_neurons(i) = sum(is_chan);
end
%
keep_mask = n_neurons > 15; % only use days when we recorded at least 15 neurons
%
unique_idx = unique_idx(keep_mask);

peak_amps_stim = zeros(numel(exp_amp_freq_data.nonstim_chan{1}.spikeTrialTimes) , numel(unique_idx)); % condition, stimulation channel
peak_amps_base = zeros(size(peak_amps_stim));
peak_amps_fr = zeros(size(peak_amps_stim));

ISI_stim = zeros(numel(exp_amp_freq_data.nonstim_chan{1}.spikeTrialTimes), numel(unique_idx), numel(ISI_edges)-1);
ISI_base = zeros(size(ISI_stim));
ISI_sim = zeros(size(ISI_stim));

for i_cond = 1:numel(exp_amp_freq_data.nonstim_chan{1}.spikeTrialTimes) % for each condition
    n_pulses = exp_amp_freq_data.nonstim_chan{1}.STIM_PARAMETERS(i_cond).nPulsesTrain;
    chan_counter = 1;
    for i_chan = unique_idx
        
        array_data_use = exp_amp_freq_data.nonstim_chan(exp_amp_freq_data.nonstim_idx == i_chan);
        % extract all spikes and trial idx
        spike_ts = [];
        spike_trial = [];
        mean_fr_train = [];
        for i_arr = 1:numel(array_data_use)
            spike_ts = [spike_ts, array_data_use{i_arr}.spikeTrialTimes{i_cond}];
            spike_trial = [spike_trial, array_data_use{i_arr}.stimData{i_cond}];
            
            % get FR during train
            train_mask = array_data_use{i_arr}.spikeTrialTimes{i_cond} > 0 & array_data_use{i_arr}.spikeTrialTimes{i_cond} < 4;
            mean_fr_train(i_arr) = sum(train_mask)/max(array_data_use{i_arr}.stimData{i_cond})/(4);

        end

        % for each trial, convolve spikes across all channels (spike_ts)
        % with a gaussian, compute peak to peak amplitudes after each pulse
        % and take mean. Also compute peak to peak amplitudes at random
        % segments during baseline (before the train)
        for i_trial = 1:max(spike_trial) % 8 trials typically
            spike_mask = spike_trial == i_trial;

            % bin spike data, then convolve with gaussian
%             [bin_spikes] = histcounts(spike_ts(spike_mask==1), bin_edges);
%             bin_spikes = bin_spikes;
%             
%             spike_train = conv(kernel,bin_spikes)./nm'/numel(array_data_use);
%             spike_train = spike_train(kernel_hl+1:end-kernel_hl);

            % compute peak to peak amplitude in window after each pulse,
            % add all up and divide by total after for loop. 
            % do same for random windows during baseline
% 
%             for i_pulse = 1:n_pulses
%                 % get pulse time
%                 pulse_time = array_data_use{1}.PULSE_TIMES{i_cond}{i_trial}(i_pulse);
%                 % get max in window after pulse_time
%                 window_idx = [find(bin_edges > pulse_time+window(1),1,'first'), find(bin_edges > pulse_time+window(2),1,'first')];
%                 peak_amps_stim(i_cond,chan_counter) = peak_amps_stim(i_cond,chan_counter) + max(spike_train(window_idx(1):window_idx(2)));
% 
%                 % pick random baseline time (-2, -0.2) s before stim
%                 rand_time = rand()*1.8 - 2;
%                 window_idx = [find(bin_edges > rand_time+window(1),1,'first'), find(bin_edges > rand_time+window(2),1,'first')];
%                 peak_amps_base(i_cond,chan_counter) = peak_amps_base(i_cond,chan_counter) + max(spike_train(window_idx(1):window_idx(2)));
%                 
%             end
            % for each trial, generate ISI histogram during baseline and during
            % the stim pulse
            spike_ts_sort = sort(spike_ts(spike_mask==1)); % use spikes from this trial
            spike_diff = diff(spike_ts_sort);
            dur_stim_mask = spike_ts_sort(1:end-1) > 0 & spike_ts_sort(1:end-1) < 4;
            before_stim_mask = spike_ts_sort(1:end-1) < 0;
            
            [ISI_temp] = histcounts(spike_diff(dur_stim_mask==1), ISI_edges);
            ISI_stim(i_cond,chan_counter,:) = ISI_stim(i_cond,chan_counter,:) + reshape(ISI_temp,1,1,size(ISI_stim,3));
            
            [ISI_temp] = histcounts(spike_diff(before_stim_mask==1), ISI_edges);
            ISI_base(i_cond,chan_counter,:) = ISI_base(i_cond,chan_counter,:) + reshape(ISI_temp,1,1,size(ISI_base,3));
            
        end
        peak_amps_stim(i_cond,chan_counter) = peak_amps_stim(i_cond,chan_counter)/(max(spike_trial)); % divide by number of pulses total
        peak_amps_base(i_cond,chan_counter) = peak_amps_base(i_cond,chan_counter)/(max(spike_trial));
    
        
        
        % simulate (poisson process) neuron's firing at rate during
        % train, compute metric for that spike train as a second
        % control
        sim_edges = 0:dt:4;
        sim_centers = sim_edges(1:end-1) + (sim_edges(2)-sim_edges(1))/2;
        for i_trial = 1:max(spike_trial)
            sim_spike_ts = [];
            for i_arr = 1:numel(mean_fr_train)
                temp_bin_spikes = poissrnd(mean_fr_train(i_arr)*dt, size(sim_centers,1), size(sim_centers,2));
                sim_spike_ts = [sim_spike_ts, sim_centers(find(temp_bin_spikes >= 1))];
            end

            sim_spike_ts = sort(sim_spike_ts);
            sim_diff = diff(sim_spike_ts);
            [ISI_temp] = histcounts(sim_diff, ISI_edges);
            ISI_sim(i_cond,chan_counter,:) = ISI_sim(i_cond,chan_counter,:) + reshape(ISI_temp,1,1,size(ISI_sim,3));
        end
%         sim_bin_spikes = sim_bin_spikes/sum(sim_bin_spikes);
        
%         sim_train = conv(kernel,sim_bin_spikes)./nm'/numel(array_data_use);
%         sim_train = sim_train(kernel_hl+1:end-kernel_hl);
%         % get max after each artificial pulse
%         for i_pulse = 1:n_pulses
%             % get pulse time
%             pulse_time = array_data_use{1}.PULSE_TIMES{i_cond}{i_trial}(i_pulse);
%             % get max in window after pulse_time
%             window_idx = [find(bin_edges > pulse_time+window(1),1,'first'), find(bin_edges > pulse_time+window(2),1,'first')];
%             peak_amps_fr(i_cond,chan_counter) = peak_amps_fr(i_cond,chan_counter) + max(sim_train(window_idx(1):window_idx(2)));
%         end
        
        peak_amps_fr(i_cond,chan_counter) = peak_amps_fr(i_cond,chan_counter)/n_pulses;
        
        % update counter
        chan_counter = chan_counter + 1;
    end

end



%% visualize ISI data
    % get ISI's < 0.5 ms
    max_ISI_edge = find(ISI_edges <= 0.5/1000, 1,'last');
    n_ISI_stim = sum(ISI_stim(:,:,1:max_ISI_edge),3);
    n_ISI_base = sum(ISI_base(:,:,1:max_ISI_edge),3);
    n_ISI_sim = sum(ISI_sim(:,:,1:max_ISI_edge),3);
        
    % get mean # ISI's < 0.5 ms, and std.
    mean_ISI_stim = mean(n_ISI_stim, 2);
    std_ISI_stim = std(n_ISI_stim, 0, 2);
    
    mean_ISI_base = mean(n_ISI_base,2);
    std_ISI_base = std(n_ISI_base, 0, 2);
    
    mean_ISI_sim = mean(n_ISI_sim,2);
    std_ISI_sim = std(n_ISI_sim,0,2);
    
    f=figure(); hold on;
    f.Position = [600 433 427 420];
    offset = [-2.5,-1.5,-0.5]*1.25;
    control_offset = [0.5,1.5,2.5]*1.25;
    for i_amp = 1:3
        idx = find(amp_data == amp_data(i_amp));
        % plot actual data
        errorbar(freq_data(idx)+offset(i_amp), mean_ISI_stim(idx),std_ISI_stim(idx),...
            'linestyle','none','marker','.','markersize',24,'color',colors(i_amp,:),'linewidth',1.2);
        
        % plot simulated control
        errorbar(freq_data(idx)+control_offset(i_amp), mean_ISI_sim(idx),std_ISI_sim(idx),...
            'linestyle','none','marker','s','markersize',6,'color',colors(i_amp,:),'linewidth',1.2);
        
    end
    
    xlabel('Frequency (Hz)');
    ylabel('Number of interspike intervals < 0.5 ms')

    xlim([40,141])
%     ylim([0,0.035])
    ax = gca;
    ax.XTick = [51,80,104,131];
    formatForLee(gcf);
    set(gca,'fontsize',14);
    ax.XMinorTick = 'off'; 
    
% statistics: compare effects of amp freq and simulated vs cortical
    amp_data = [20,40,60,20,40,60,20,40,60,20,40,60];
    freq_data = [131,131,131,104,104,104,80,80,80,51,51,51];
    mean_n_stim = mean(n_ISI_stim,2);    
    mean_n_sim = mean(n_ISI_base,2);
    
    amp_all = [];
    freq_all = [];
    is_cortical = [];
    peak_all = [];
    
    for i = 1:size(n_ISI_stim,1)
        for j = 1:size(n_ISI_stim,2)
            % add stim data
            amp_all(end+1) = amp_data(i);
            freq_all(end+1) = freq_data(i);
            is_cortical(end+1) = 1;
            peak_all(end+1) = n_ISI_stim(i,j);

            % add simulated data
            amp_all(end+1) = amp_data(i);
            freq_all(end+1) = freq_data(i);
            is_cortical(end+1) = 0;
            peak_all(end+1) = n_ISI_sim(i,j);
        end
    end

    % build table
    data_table = table(peak_all',categorical(is_cortical)',amp_all',freq_all','VariableNames',...
        {'peak','is_cort','amp','freq'});
    mdlspec = 'peak ~ is_cort + amp*freq';
    
    mdl = fitlm(data_table,mdlspec)
    
%% make raster plot of 1 trial, each row is a different neuron

% load map files
    map_data_han = loadMapFile('R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp');
    map_data_dunc = loadMapFile('R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp');
    
    idx = 5; % or 8 should have the most neurons
    i_cond = 3; % 60 uA, 131 Hz
    i_trial = 2; 
    
    array_data_use = exp_amp_freq_data.nonstim_chan(exp_amp_freq_data.nonstim_idx == unique_idx(idx));
    
    % plot neurons in order of distance from stim chan
    dists_from_stim = [];
    for i_arr = 1:numel(array_data_use)
        % pick map data
        if(strcmpi(array_data_use{i_arr}.monkey, 'Han')==1)
            map_data = map_data_han;
        else
            map_data = map_data_dunc;
        end
        
        stim_chan_idx = find(map_data.chan == array_data_use{i_arr}.CHAN_SENT{1});
        rec_chan_idx = find(map_data.chan == array_data_use{i_arr}.CHAN_REC);
        
        stim_chan_pos = [map_data.row(stim_chan_idx), map_data.col(stim_chan_idx)];
        rec_chan_pos = [map_data.row(rec_chan_idx), map_data.col(rec_chan_idx)];
        
        dists_from_stim(i_arr) = sqrt(sum((stim_chan_pos-rec_chan_pos).^2));
    end
    
    [~,row_idx] = sort(dists_from_stim);
    
    % extract all spikes from a given trial, set trial idx to i_arr idx
    spike_ts = [];
    spike_trial = [];
    
    for i_arr = 1:numel(array_data_use)
        spike_mask = array_data_use{i_arr}.stimData{i_cond} == i_trial;
        spike_ts = [spike_ts, array_data_use{i_arr}.spikeTrialTimes{i_cond}(spike_mask)];
        spike_trial = [spike_trial, zeros(1,sum(spike_mask))+i_arr];
    end

    % plot raster, each row is a different neuron (stored in spike trial)
    opts_plot = []; opts_save = [];
    opts_plot.X_LABEL = 'Time after train onset (s)';
    opts_plot.Y_LABEL = '';
    opts_plot.X_LIMITS = [-100,200]/1000; % in s
    opts_plot.Y_LIMITS = [0.5,numel(array_data_use)+0.5];
    opts_plot.MARKER_STYLE = 'line';
    
    plotRasterLIB(spike_ts, spike_trial, opts_plot, opts_save);
    hold on;
    f=gcf;
    f.Position = [600 433 462 420];
    % plot stim times
    stim_times = array_data_use{1}.PULSE_TIMES{i_cond}{i_trial};

    for i = 1:numel(stim_times)
        plot([stim_times(i), stim_times(i)],opts_plot.Y_LIMITS,'r-','linewidth',0.25)
    end
    
%% plot example ISI

    i_chan = 1;
    i_cond = 3; % 60uA, 131 Hz
    
    f=figure('Position',[600 433 420 420]); hold on;    
    ISI_centers = 1000*(ISI_edges(1:end-1) + mode(diff(ISI_edges))/2);
    plot(ISI_centers, squeeze(ISI_sim(i_cond,i_chan,:)),'s','color',getColorFromList(1,0),'markersize',10);
    
    plot(ISI_centers,squeeze(ISI_stim(i_cond,i_chan,:)),'.','color',getColorFromList(1,1),'markersize',24);
    
    xlabel('Time between spikes across array (ms)');
    ylabel('Number of occurences');
    
    formatForLee(gcf);
    set(gca,'fontsize',14);
    xlim([0,3])

    
    
%% bin spike data, then convolve with gaussian
    [bin_spikes] = histcounts(spike_ts, bin_edges);

    spike_train = conv(kernel,histcounts(spike_ts,bin_edges))./nm'/numel(array_data_use);
    spike_train = spike_train(kernel_hl+1:end-kernel_hl);
    
    f=figure('Position',[600 433 467 170]); hold on;
    
    plot(bin_edges(1:end-1)+mode(diff(bin_edges))/2,spike_train,'k-','linewidth',2);
    xlim(opts_plot.X_LIMITS);
    for i = 1:numel(stim_times)
        plot([stim_times(i), stim_times(i)],[0,0.06],'r-','linewidth',0.25)
    end
    xlabel('Time after train onset (s)');
    
    set(gca,'fontsize',14)
    formatForLee(gcf);
    set(gca,'YTickLabel',[]);

    
%% visualize convolution data (points with errorbars)

    colors = inferno(4);

    f=figure('Position',[680 558 356 420]); hold on
    f.Name = 'nonstim_chan_sync';
    
    mean_peak_stim = mean(peak_amps_stim,2);
    std_peak_stim = std(peak_amps_stim,0,2);
    
    mean_peak_fr = mean(peak_amps_fr,2);
    std_peak_fr = std(peak_amps_fr,0,2);
    
    offset = [-2.5,-1.5,-0.5]*1.25;
    control_offset = [0.5,1.5,2.5]*1.25;
    for i_amp = 1:3
        idx = find(amp_data == amp_data(i_amp));
        % plot actual data
        errorbar(freq_data(idx)+offset(i_amp), mean_peak_stim(idx),std_peak_stim(idx),...
            'linestyle','none','marker','.','markersize',24,'color',colors(i_amp,:),'linewidth',1.2);
        
        % plot simulated control
        errorbar(freq_data(idx)+control_offset(i_amp), mean_peak_fr(idx),std_peak_fr(idx),...
            'linestyle','none','marker','s','markersize',6,'color',colors(i_amp,:),'linewidth',1.2);
        
    end
    
    xlabel('Frequency (Hz)');
    ylabel('Mean peak amplitude (a.u.)')

    xlim([40,141])
%     ylim([0,0.035])
    ax = gca;
    ax.XTick = [51,80,104,131];
    formatForLee(gcf);
    set(gca,'fontsize',14);
    ax.XMinorTick = 'off'; 
