%% get experimental data -- stim channel response
    exp_input_data.home_computer = 1;
    exp_input_data.amp_list = [5,10,15,20,25,30,40,50,100];
    
    [exp_data] = getExperimentStimChannelData(exp_input_data);
    exp_array_data = exp_data.array_data;
    exp_array_data = adjustArrayDataSpikeTimes(exp_array_data, 0.453/1000); % stim pulse length
    exp_array_data = getBaselineFiringRate(exp_array_data,[-80,-10]/1000); % window relative to stim onset
 
%% get amplitude list for each neuron
    amps = {};
    num_stims = [];
    for i_unit = 1:numel(exp_array_data)
        amps{i_unit,1} = [exp_array_data{i_unit}.STIM_PARAMETERS.amp1];
        num_stims = [num_stims, exp_array_data{i_unit}.numStims];
    end
    
    
%% plot raster of example neuron

    raster_input_data = [];
    raster_input_data.x_lim = [-10,10]; % ms
    raster_input_data.amp_list = exp_input_data.amp_list;
    raster_input_data.marker_style = 'line'; % line is the correct way, but much slower
    raster_input_data.plot_amp = 1;
    
    load('D:\Lab\Data\stim_ephys_paper\artifact_analysis\sim_study\duke_blank_times');
    raster_input_data.duke_blank_times = blank_times;
    load('D:\Lab\Data\stim_ephys_paper\artifact_analysis\sim_study\black_blank_times');
    raster_input_data.black_blank_times = blank_times;
    
    for exp_idx = 5%1:numel(exp_array_data)%5 %4, 5, 13, 25 (14,12)
        raster_input_data.is_model = 0;
        plotModelExpAmpRaster(exp_array_data{exp_idx},raster_input_data);
    end

%% Activation threshold
    activation_input_data.spike_window = [0.5,5]/1000;
    
    activation_input_data.remove_intrinsic = 1;
    activation_input_data.amp_list = exp_input_data.amp_list;
    
    activation_input_data.use_stats = 1; % or use percent cutoff. Threshold is used for both
    
    if(activation_input_data.use_stats) 
        activation_input_data.sub_baseline = 0;
        activation_input_data.threshold = 0.05;
    else
        activation_input_data.sub_baseline = 1;
        activation_input_data.threshold = 0.5;
    end
  
    activation_input_data.test_thresholds = 0:0.05:1;
    
    activation_input_data.is_model = 0;
    exp_threshold_data = getActivationThreshold(exp_array_data,activation_input_data);

%% plot activation threshold data
    f=figure('Position',[2918,452,250,327]); hold on;
    f.Name = 'Han_Duncan_activationThreshold';
        
    for i_amp = 1:numel(exp_input_data.amp_list)
        num_neurons(i_amp) = sum(exp_threshold_data.thresholds == exp_input_data.amp_list(i_amp));
        
        bar(exp_input_data.amp_list(i_amp), num_neurons(i_amp),4.5,...
            'FaceColor',[0.2,0.2,0.2]);
        
    end
    
    formatForLee(gcf); set(gca,'fontsize',14);
    ylabel('Number of neurons');
    xlabel('Amplitude (\muA)');
    xlim([0,30]);
    set(gca,'XTick',[0:10:40],'XTickLabel',{'0','10','20','30'});
    ax = gca; ax.XAxis.MinorTickValues = [0:5:45];
    ax.YAxis.MinorTickValues = [0:1:18];
    
%% plot proportion responsive across amplitudes
    f=figure('Position',[2255 443 450 420]); hold on;
    f.Name = 'Han_Duncan_singlePulse_proportionResponsive';
    amp_list = [5,10,15,20,25,30,40,50,100];
    boxplot_params = [];
    boxplot_params.use_same_color_for_all = 1; % omits median color
    
    boxplot_params.median_color = 'k';
    boxplot_params.box_width = 3;
    boxplot_params.linewidth = 1.5;
    
    % variables for stats at bottom of section
    amp_tbl = [];
    id_tbl = [];
    num_spikes_tbl = [];
    for i_amp = 1:numel(exp_input_data.amp_list)        
        boxplot_wrapper(amp_list(i_amp),exp_threshold_data.avg_num_spikes(:,i_amp),boxplot_params);
        
        num_add = sum(~isnan(exp_threshold_data.avg_num_spikes(:,i_amp)));
        amp_tbl = [amp_tbl; exp_input_data.amp_list(i_amp)*ones(num_add,1)];
        id_tbl = [id_tbl; find(~isnan(exp_threshold_data.avg_num_spikes(:,i_amp)))];
        num_spikes_tbl = [num_spikes_tbl; exp_threshold_data.avg_num_spikes(~isnan(exp_threshold_data.avg_num_spikes(:,i_amp)),i_amp)];
    end
    
    % 0uA case
%     boxplot_wrapper(0,exp_threshold_data.avg_num_spikes_null,boxplot_params);
    
    formatForLee(gcf); set(gca,'fontsize',14);
    ylabel('Spikes per pulse above baseline');
    xlabel('Amplitude (\muA)');
    xlim([0,105]);
    
    set(gca,'XTick',[0:10:100]);
    ax = gca; ax.XAxis.MinorTickValues = [0:5:100];
    ax.YAxis.MinorTickValues = [0:1:18];

    
    
    num_tbl = table(amp_tbl, id_tbl,num_spikes_tbl,...
        'VariableNames',{'amp','neuron','num_spikes'});
    num_tbl.neuron = categorical(num_tbl.neuron);

    num_tbl_spec = 'num_spikes~amp+neuron';
    num_mdl = fitlm(num_tbl,num_tbl_spec)

%% magnitude and latency of volleys
    lat_input_data.amp_list = exp_input_data.amp_list;
    lat_input_data.peak_window = [-10,10]/1000; % s
    lat_input_data.bin_size = 0.1/1000; % s
    
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
    
    f=figure('Position',[2098,350,350,420]);
    f.Name = 'Han_Duncan_changePeakLatency';
    boxplot_params.master_color = [0.4,0.4,0.4];
    x_data = []; y_data = [];
    
    delta_amp = exp_latency_data.delta_amp - exp_latency_data.delta_min_amp;
    unique_amp_list = unique(delta_amp);
    for i_amp = 1:numel(unique_amp_list)
        data_mask = delta_amp == unique_amp_list(i_amp);
        x_data = [x_data;unique_amp_list(i_amp)+zeros(sum(data_mask),1)];
        y_data = [y_data;exp_latency_data.delta_lat(data_mask)'*1000];
        if(sum(data_mask) > 2)
            boxplot_wrapper(unique_amp_list(i_amp),exp_latency_data.delta_lat(data_mask)*1000,boxplot_params);
%             hold on
%             plot(unique_amp_list(i_amp)+2*rand(sum(data_mask),1)-1,exp_latency_data.delta_lat(data_mask)*1000,'ko','markersize',6);
        end
    end
    
    
    xlabel('\Delta Amplitude (\muA)');
    ylabel('\Delta peak latency (ms)');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    xlim([0,45])
    
    lat_mdl = fitlm(x_data,y_data,'Intercept',false)
%     lat_mdl = LinearModel.fit(x_data-min(x_data),y_data,'Intercept',false)
    
%% latency vs. standard dev of spikes in a peak

    f=figure('Position',[2885 476 300 327]); hold on
    f.Name = 'Han_duncan_stimChan_spikeTimingVsSpread';
    % latency vs. width
    amps_plot = [2,5,8,9];
    color_list = inferno(numel(amps_plot)+1); % remove most yellow color
    unique_amps = unique(exp_latency_data.amp_idx);
    color_cnt=1;
    lat_list = []; std_list = []; amp_list = []; id_list = []; group_id_list = [];
    for i_amp = 1:numel(unique_amps)
            mask = exp_latency_data.amp_idx == unique_amps(i_amp);
            lat_list = [lat_list; exp_latency_data.lat_list(mask)'*1000];
            std_list = [std_list; exp_latency_data.std_list(mask)'*1000];
            amp_list = [amp_list; unique_amps(i_amp)*ones(sum(mask),1)];
            id_list = [id_list; exp_latency_data.unit_idx(mask)'];
            group_id_list = [group_id_list; exp_latency_data.peak_num(mask)'];
            
        if(any(i_amp==amps_plot))
            plot((exp_latency_data.lat_list(mask)*1000),(exp_latency_data.std_list(mask)*1000),...
                'o','color',color_list(color_cnt,:),'markersize',4)
            color_cnt = color_cnt+1;
        end
    end
    xlabel('Peak latency (ms)');
    ylabel('Spike timing std (ms)');
    formatForLee(gcf);
    set(gca,'fontsize',14)
    
    mask = amp_list <= Inf;
    std_tbl = table(lat_list(mask), std_list(mask), id_list(mask), amp_list(mask), group_id_list(mask),...
        'VariableNames',{'lat','std','neuron','amp','group_num'});
    std_tbl.neuron = categorical(std_tbl.neuron);
    
    std_mdl_spec = 'std~lat+amp+neuron+group_num';
    std_mdl = fitlm(std_tbl,std_mdl_spec)
    
    lat_mdl_spec = 'lat~amp:group_num + group_num +neuron';
    lat_mdl = fitlm(std_tbl,lat_mdl_spec)
    
    
    
%% num peaks vs amplitude

    f=figure('Position',[740 476 678 327]); hold on;
    f.Name = 'Han_duncan_numPeaksVsAmplitude';

    offset = [-1.25,-0.75,-0.25,0.25,0.75,1.25]*2;
    red_fact = 0.75;
    
    num_peak_data = zeros(numel(lat_input_data.amp_list),1+max([exp_latency_data.num_peaks_amp(:,1);exp_latency_data.num_peaks_amp_assume(:,1)]));
    percent_neurons = [];
    for i_amp = 1:numel(lat_input_data.amp_list)
        data_mask = exp_latency_data.num_peaks_amp(:,2) == lat_input_data.amp_list(i_amp);
        
        for i_num = 1:size(num_peak_data,2)
            num_peaks = i_num-1;
            num_peak_data(i_amp,i_num) = sum(exp_latency_data.num_peaks_amp(data_mask,1) == num_peaks);

            percent_neurons = num_peak_data(i_amp,i_num)/sum(data_mask);
            
%             % draw square for this peak num and amplitude
            color_use = [1,1 - percent_neurons*red_fact,1 - percent_neurons*red_fact]; % red only map
            rectangle('Position',[lat_input_data.amp_list(i_amp)-2.5,num_peaks-0.5,5,1],'EdgeColor','k','FaceColor',color_use)
            % write number of units with that number of peaks in square
            text(lat_input_data.amp_list(i_amp),num_peaks,num2str(round(percent_neurons*100,0)),'horizontalalignment','center','fontsize',12)
        end
    end
    colormap([ones(100,1,1), ones(100,1,1)-linspace(0,red_fact,100)', ones(100,1,1)-linspace(0,red_fact,100)']);

    b=colorbar;
    b.Label.String = 'Percent of neurons';
    
    % plot num_peaks_assume for 100uA
    for i_num = 1:size(num_peak_data,2)
        num_peaks = i_num-1;
        num_neurons = size(exp_latency_data.num_peaks_amp_assume,1);
        num_peak_data(1,i_num) = sum(exp_latency_data.num_peaks_amp_assume(:,1) == num_peaks);

        percent_neurons = num_peak_data(1,i_num)/num_neurons;
        
        % draw square for this peak num and amplitude
        color_use = [1-percent_neurons*red_fact,1 - percent_neurons*red_fact,1 - percent_neurons*red_fact]; % red only map
        rectangle('Position',[5+lat_input_data.amp_list(end)-2.5,num_peaks-0.5,5,1],'EdgeColor','k','FaceColor',color_use)
        % write number of units with that number of peaks in square
        text(lat_input_data.amp_list(end)+5,num_peaks,num2str(round(percent_neurons*100,0)),'horizontalalignment','center','fontsize',12)
    end
    
    xlabel('Amplitude (\muA)');
    ylabel('Num peaks');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    
    % set y and x tick marks
    ax=gca;
    ax.YTick = 0:2:max(exp_latency_data.num_peaks_amp(:,1));
    ax.YRuler.MinorTickValues = 0:1:max(exp_latency_data.num_peaks_amp(:,1))+1;
    ax.XTick = [5,15,25,40,50,100];
    ax.XRuler.MinorTickValues = lat_input_data.amp_list;
    
    ylim([-0.5,size(num_peak_data,2)]);
    xlim([0,110])
    
    % setup linear model for statistics. then build linear model
%     amp_mask = exp_latency_data.num_peaks_amp(:,2) ~= 100;
    peak_tbl = table(exp_latency_data.num_peaks_amp(:,2),exp_latency_data.num_peaks_amp(:,1),exp_latency_data.num_peaks_id(:),...
        'VariableNames',{'amp','peaks','id'});
    peak_tbl.id = categorical(peak_tbl.id);
    
    mdl_spec = 'peaks~amp+id';
    peak_mdl = fitlm(peak_tbl,mdl_spec)
    
%% compute inhibition duration and plot across amplitudes

    inhib_input_data = [];
    inhib_input_data.pre_window = [-80,-10]/1000; % s 
    inhib_input_data.post_window = [0,220]/1000; % s
    inhib_input_data.bin_window = [inhib_input_data.pre_window(1),inhib_input_data.post_window(2)+10/1000];
    inhib_input_data.max_time_start = 40/1000; % s
    inhib_input_data.bin_size = 5/1000; % s
    inhib_input_data.kernel_length = 1;
    inhib_input_data.blank_time = [0,10]/1000; % s
    
    inhib_input_data.num_consec_bins = 2;
    inhib_input_data.amp_list = [5,10,15,20,25,30,40,50,100];
    
    exp_array_data = rebinArrayData(exp_array_data,inhib_input_data.bin_size*1000);
    exp_inhib_data = getInhibitionDurationAmpWrapper(exp_array_data,inhib_input_data);
    
%% plot percent of cells that have inhibition, and duration at each
    % amplitude
    boxplot_params = [];
    boxplot_params.use_same_color_for_all = 1;
    boxplot_params.master_color = [0.4,0.4,0.4];
    boxplot_params.median_color = 'k';
    boxplot_params.box_width = 2.4;
    boxplot_params.linewidth = 2;
    
    f=figure('Position',[2294 486 737 420]);
    f.Name = 'Han_duncan_inhibition';
    subplot(1,2,1)
    offset = [-1,1];
    amp_data = []; inhib_dur_data = []; id_list = [];
    for i_amp = 1:numel(inhib_input_data.amp_list)   
        data = exp_inhib_data.inhib_dur(:,i_amp)*1000;
        temp_id = 1:1:numel(data)';
        temp_id(isnan(data)) = [];
        data(isnan(data)) = [];
        
        amp_data = [amp_data; inhib_input_data.amp_list(i_amp)*ones(numel(data),1)];
        inhib_dur_data = [inhib_dur_data; data];
        id_list = [id_list; temp_id'];
        boxplot_wrapper(inhib_input_data.amp_list(i_amp), data, boxplot_params);

    end
    % format
    formatForLee(gcf)
    set(gca,'fontsize',14)
    ylabel('Inhibition duration (ms)');
    xlabel('Amp (\muA)');
    xlim([0,110]);
    
    % plot percent inhib
    subplot(1,2,2); hold on; % percent inhib
    percent_inhib = sum(~isnan(exp_inhib_data.inhib_dur),1)/size(exp_inhib_data.inhib_dur,1);
    plot(inhib_input_data.amp_list,percent_inhib,'color','k','linewidth',1.5,'marker','.','markersize',14)

    % format
    formatForLee(gcf)
    set(gca,'fontsize',14)
    ylabel('Fraction of cells with inhibitory response');
    xlabel('Amp (\muA)');
    ylim([0,1])
    xlim([0,110]);
    
    
    % inhibition duration stats
    inhib_tbl = table(amp_data,inhib_dur_data,categorical(id_list),'VariableNames',{'amp','dur','id'});
    
    mdlspec = 'dur~amp + id';
    inhib_mdl = fitlm(inhib_tbl,mdlspec)
    
    
    
    