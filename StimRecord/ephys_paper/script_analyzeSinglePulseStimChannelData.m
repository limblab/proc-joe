%% get experimental data -- stim channel response
    exp_input_data.home_computer = 1;
    exp_input_data.amp_list = [5,10,15,20,25,30,40,50,100];
    
    [exp_data] = getExperimentStimChannelData(exp_input_data);
    exp_array_data = exp_data.array_data;
    exp_array_data = adjustArrayDataSpikeTimes(exp_array_data, 0.453/1000); % stim pulse length
    exp_array_data = getBaselineFiringRate(exp_array_data,[-80,-5]/1000); % window relative to stim onset
 
    
%% plot raster of example neuron

    raster_input_data = [];
    raster_input_data.x_lim = [-15,30]; % ms
    raster_input_data.amp_list = exp_input_data.amp_list;
    raster_input_data.marker_style = '.'; % line is the correct way, but much slower
    raster_input_data.plot_amp = 1;
    
    
    for exp_idx = 1:numel(exp_array_data)%5 %4, 5, 13, 25
        raster_input_data.is_model = 0;
        plotModelExpAmpRaster(exp_array_data{exp_idx},raster_input_data);
    end

%% Activation threshold
    activation_input_data.spike_window = [1,5]/1000;
    activation_input_data.remove_intrinsic = 1;
    activation_input_data.sub_baseline = 1;
    activation_input_data.amp_list = exp_input_data.amp_list;
    activation_input_data.threshold = 0.5;
    activation_input_data.test_thresholds = 0:0.05:1;
    
    activation_input_data.is_model = 0;
    exp_threshold_data = getActivationThreshold(exp_array_data,activation_input_data);

% plot activation threshold data
    f=figure('Position',[2918,452,446,420]);
    f.Name = 'Han_Duncan_activationThreshold';
    subplot(1,2,1)
    % experimental data
    boxplot_params = [];
    boxplot_params.use_same_color_for_all = 1;
    boxplot_params.master_color = 'k';
%     mask = exp_threshold_data.is_responsive == 1;
%     data_all(end+1:end+sum(mask)) = exp_threshold_data.thresholds(mask);
%     group_all(end+1:end+sum(mask)) = 4;
%     
%     vs = violinplot(data_all,group_all)
    data = exp_threshold_data.thresholds(exp_threshold_data.is_responsive==1);
    boxplot_wrapper(4, data, boxplot_params);
    
    formatForLee(gcf); set(gca,'fontsize',14);
    ylabel('Activation threshold (\muA)');
    set(gca,'XTick',[])
%
    % plot percent responsive all cells against different thresholds (cdf)
    
    percent_resp_thresh = zeros(size(activation_input_data.test_thresholds));
    for i_thresh = 1:numel(activation_input_data.test_thresholds)
        percent_resp_thresh(i_thresh) = sum(max(exp_threshold_data.percent_responsive,[],2) > activation_input_data.test_thresholds(i_thresh));
    end
    percent_resp_thresh = percent_resp_thresh./size(exp_threshold_data.percent_responsive,1);
    
    subplot(1,2,2); hold on
    % experimental data
%     perc_resp = sum(exp_threshold_data.is_responsive)/numel(exp_threshold_data.is_responsive);
%     b=bar(4, perc_resp);
%     b.FaceColor = 'k';
    plot(activation_input_data.test_thresholds, percent_resp_thresh,'k','linewidth',2);
    plot([0.5,0.5],[0,percent_resp_thresh(find(activation_input_data.test_thresholds==0.5))],'k--','linewidth',1)
    plot([0,0.5],[0,0]+percent_resp_thresh(find(activation_input_data.test_thresholds==0.5)),'k--','linewidth',1)
    formatForLee(gcf); set(gca,'fontsize',14);
    ylabel('% cells responsive above threshold');
    xlabel('Threshold (% pulses)')


%% magnitude and latency of volleys
    lat_input_data.amp_list = exp_input_data.amp_list;
    lat_input_data.peak_window = [0,10]/1000; % s
    lat_input_data.bin_size = 0.2/1000; % s
    
    lat_input_data.use_gauss_filter = 1;
    lat_input_data.dt = 0.0001;
    lat_input_data.kernel_SD = 0.0002;
    
    lat_input_data.is_model = 0;
    exp_latency_data = getLatencyOfPeaks(exp_array_data,lat_input_data);

    
%% plot change in latency of first peak vs. amplitude (peak timing across amps)
    % change in latency relative to peak at lowest amplitude for each
    % unit
    
    boxplot_params = [];
    boxplot_params.use_same_color_for_all = 1; % omits median color
    
    boxplot_params.median_color = 'k';
    boxplot_params.box_width = 3;
    
    f=figure('Position',[2098,350,429,420]);
    f.Name = 'Han_Duncan_changePeakLatency';
    boxplot_params.master_color = 'k';

    for i_amp = 1:numel(lat_input_data.amp_list)
        data_mask = exp_latency_data.delta_amp == lat_input_data.amp_list(i_amp);

        if(sum(data_mask) > 2)
            boxplot_wrapper(lat_input_data.amp_list(i_amp),exp_latency_data.delta_lat(data_mask)*1000,boxplot_params);
        end
    end
    
    xlabel('Amplitude (\muA)');
    ylabel('\Delta peak latency (ms)');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    xlim([0,110])
    
%% latency vs. standard dev of spikes in a peak

    f=figure();
    f.Name = 'Han_duncan_stimChan_spikeTimingVsSpread'
    % latency (or peak num) vs. width
    plot(exp_latency_data.lat_list*1000,exp_latency_data.std_list*1000,'k.')
    xlabel('Peak latency (ms)');
    ylabel('Spike timing std (ms)');
    formatForLee(gcf);
    set(gca,'fontsize',14)
    
%% num peaks vs amplitude


    f=figure('Position',[740 476 678 247]); hold on;
    f.Name = 'Han_duncan_numPeaksVsAmplitude';

    red_fact = 0.75;
    
    num_peak_data = zeros(numel(lat_input_data.amp_list),1+max(exp_latency_data.num_peaks_amp(:,1)));
    for i_amp = 1:numel(lat_input_data.amp_list)
        data_mask = exp_latency_data.num_peaks_amp(:,2) == lat_input_data.amp_list(i_amp);
        
        for i_num = 1:size(num_peak_data,2)
            num_peaks = i_num-1;
            num_peak_data(i_amp,i_num) = sum(exp_latency_data.num_peaks_amp(data_mask,1) == num_peaks);

            % draw square for this peak num and amplitude
            color_use = [1,1 - num_peak_data(i_amp,i_num)*red_fact/sum(data_mask),1 - num_peak_data(i_amp,i_num)*red_fact/sum(data_mask)]; % red only map
            rectangle('Position',[lat_input_data.amp_list(i_amp)-2.5,num_peaks-0.5,5,1],'EdgeColor','k','FaceColor',color_use)
            % write number of units with that number of peaks in square
            text(lat_input_data.amp_list(i_amp),num_peaks,num2str(num_peak_data(i_amp,i_num)),'horizontalalignment','center','fontsize',12)
        end
    end
    
    xlabel('Amplitude (\muA)');
    ylabel('Num peaks');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    
    % set y and x tick marks
    ax=gca;
    ax.YTick = 0:2:max(exp_latency_data.num_peaks_amp(:,1));
    ax.YRuler.MinorTickValues = 0:1:max(exp_latency_data.num_peaks_amp(:,1));
    ax.XTick = [5,15,25,40,50,100];
    ax.XRuler.MinorTickValues = lat_input_data.amp_list;
    
    ylim([-0.5,max(exp_latency_data.num_peaks_amp(:,1))+0.5]);
    xlim([0,110])
    
%% compute inhibition duration and plot across amplitudes

    inhib_input_data = [];
    inhib_input_data.pre_window = [-80,-5]/1000; % s 
    inhib_input_data.post_window = [0,220]/1000; % s
    inhib_input_data.bin_window = [inhib_input_data.pre_window(1),inhib_input_data.post_window(2)+10/1000];
    inhib_input_data.max_time_start = 40/1000; % s
    inhib_input_data.bin_size = 5/1000; % s
    inhib_input_data.kernel_length = 2;
    inhib_input_data.blank_time = [0,10]/1000; % s
    
    inhib_input_data.num_consec_bins = 2;
    inhib_input_data.amp_list = [5,10,15,20,25,30,40,50,100];
    
    exp_array_data = rebinArrayData(exp_array_data,inhib_input_data.bin_size*1000);
    exp_inhib_data = getInhibitionDurationAmpWrapper(exp_array_data,inhib_input_data);
    
%% plot percent of cells that have inhibition, and duration at each
    % amplitude
    boxplot_params = [];
    boxplot_params.use_same_color_for_all = 1;
    boxplot_params.master_color = 'k';
    boxplot_params.box_width = 2.4;
    boxplot_params.linewidth = 2;
    
    f=figure('Position',[2294 486 737 420]);
    f.Name = 'Han_duncan_inhibition';
    subplot(1,2,1)
    offset = [-1,1];
    for i_amp = 1:numel(inhib_input_data.amp_list)   
        data = exp_inhib_data.inhib_dur(:,i_amp)*1000;
        data(isnan(data)) = [];

        boxplot_wrapper(inhib_input_data.amp_list(i_amp)+offset(i_type), data, boxplot_params);

    end
    % format
    formatForLee(gcf)
    set(gca,'fontsize',14)
    ylabel('Inhib dur (ms)');
    xlabel('Amp (\muA)');
    xlim([0,110]);
    
    % plot percent inhib
    subplot(1,2,2); hold on; % percent inhib
    percent_inhib = sum(~isnan(exp_inhib_data.inhib_dur),1)/size(exp_inhib_data.inhib_dur,1);
    plot(inhib_input_data.amp_list,percent_inhib,'color','k','linewidth',1.5,'marker','.','markersize',14)

    % format
    formatForLee(gcf)
    set(gca,'fontsize',14)
    ylabel('Fraction cells with inhib response');
    xlabel('Amp (\muA)');
    ylim([0,1])
    xlim([0,110]);