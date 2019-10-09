%% look at single pulse data. Describe activation thresholds, latency of
% different volleys, inhibition strength as a function of amplitude and
% other parameters.
% Need to create a cell array called arrayData where each entry is a recorded
% neuron. analyzeStimData does this


%% load in arrayData

%% switch arrayData over to stimulation offset if desired

%% rebin data if desired. Also switch to 'time after stimulation offset' if desired
    inputData.bin_size = 0.2; % in ms

    arrayData = rebinArrayData(arrayData,inputData);

%% plot rasters and psth for each condition and neuron

%     for arrIdx = 1:numel(arrayData)
    arrIdx = 1;
        % plot raster, and PSTH for the given unit above
        
    %     optsPlotFunc.BIN_SIZE = optsExtract.BIN_SIZE;
        optsPlotFunc.BIN_SIZE = mode(diff(arrayData{arrIdx}.binEdges{1,1}));
        optsPlotFunc.FIGURE_SAVE = 0;
        optsPlotFunc.FIGURE_DIR = [];
        optsPlotFunc.FIGURE_PREFIX = [arrayData{arrIdx}.monkey,'_long_'];

        optsPlotFunc.PRE_TIME = 15/1000;
        optsPlotFunc.POST_TIME = 20/1000;
        optsPlotFunc.SORT_DATA = '';

        optsPlotFunc.MARKER_STYLE  = '.';
        
        optsPlotFunc.PLOT_AFTER_STIMULATION_END = 1;
        optsPlotFunc.STIMULATION_LENGTH = [];

        rasterPlots = plotRasterStim(arrayData{arrIdx},1,optsPlotFunc);

        optsPlotFunc.PLOT_ALL_ONE_FIGURE = 0;
        optsPlotFunc.PLOT_LINE = 1;
        optsPlotFunc.PLOT_TITLE = 0;    
        optsPlotFunc.PLOT_ALL_WAVES_ONE_FIGURE = 0;
        
    %   
%         PSTHPlots = plotPSTHStim(arrayData{arrIdx},arrayData{arrIdx}.NN,optsPlotFunc);

%     end


%% get number of spikes (probability) as a function of amplitude and plot
% that for each neuron (on one plot probably)

    optsSpikesPlot.PRE_WINDOW = [-80,-5]/1000; % in s
    optsSpikesPlot.POST_WINDOW = [0,20]/1000; % in s
    optsSpikesPlot.PW1 = 200;
    optsSpikesPlot.PW2 = 200;
    optsSpikesPlot.POL = 0; % 0 is cathodic first
    optsSpikesPlot.AMP = [];
    optsSpikesPlot.PLOT_PULSE_WIDTH = 0;
    
    optsSpikesPlot.ADJUST_SPIKE_TIMES = 1;
    
    spikesStruct = {}; amp = []; spikes_evoked = []; unit_idx = [];
    
    
    for unit = 1:numel(arrayData)
        if(unit == 1)
            optsSpikesPlot.MAKE_FIGURE = 1;
        else
            optsSpikesPlot.MAKE_FIGURE = 0;
        end
        [spikesStruct{unit},figure_handles] = plotSpikesPerCondition(arrayData{unit},optsSpikesPlot);
        hold on
        
        formatForLee(gcf)
        xlabel('Amplitude (\muA)');
        ylabel('Spikes per stimulation');
        set(gca,'fontsize',14)
        
        amp = [amp,spikesStruct{unit}.amp(spikesStruct{unit}.keep_mask==1)];
        spikes_evoked = [spikes_evoked, spikesStruct{unit}.prob_spike(spikesStruct{unit}.keep_mask==1)];
        unit_idx = [unit_idx, unit + zeros(1,sum(spikesStruct{unit}.keep_mask==1))];
    end  
    
    amps_plot = unique(amp);
    mean_spikes_evoked = zeros(1,numel(amps_plot));
    
    for a = 1:numel(amps_plot)
        amp_mask = amp == amps_plot(a);
        mean_spikes_evoked(a) = mean(spikes_evoked(amp_mask));
        
    end
    
    tbl = table(amp',spikes_evoked',unit_idx','VariableNames',{'amp','spikes_evoked','unit'});
    mdl = fitlm(tbl,'spikes_evoked~amp+unit')
    
    
%% make spike time distribution plot
    amps_plot = [5,10,15,20,25,30,40,50,100];
    spike_times_amp = cell(numel(amps_plot),1); % preallocate space  
    bin_edges = [0:1:20];
    colors = inferno(numel(amps_plot) + 3);
    f_latency=figure(); hold on;
    
    for a = 1:numel(amps_plot)
        spike_times_amp{a} = [];
        for unit = 1:numel(spikesStruct)
            amp_idx = find(spikesStruct{unit}.amp == amps_plot(a) & spikesStruct{unit}.keep_mask == 1);
            if(~isempty(amp_idx))
                spike_times_amp{a} = [spike_times_amp{a}, spikesStruct{unit}.spike_times_post_stim{amp_idx}*1000];
            end
        end
        plot(bin_edges(1:end-1)+mode(diff(bin_edges)/2),histcounts(spike_times_amp{a},bin_edges),'color',colors(a,:),'linewidth',1.5)
        hold on
    end
    
    xlabel('Time after stimulation offset (ms)');
    ylabel('Number of spikes');
    formatForLee(gcf)
    set(gca,'fontsize',14)
    l=legend('5','10','15','20','25','30','40','50','100');
    set(l,'box','off','fontsize',14);
% make cdf plot for latencies (dependent on previous section)
    
    x_data = {};
    y_data = {};
    colors = inferno(numel(amps_plot) + 3);
    figure(); hold on
    for a = 1:numel(amps_plot)
                
        x_data{a} = sort([0;spike_times_amp{a}';spike_times_amp{a}';max(bin_edges)]);
        y_data_list = (0:numel(spike_times_amp{a}));
        y_data{a} = reshape(repmat(y_data_list,2,1),[1,numel(y_data_list)*2]);
        
        plot(x_data{a},y_data{a},'color',colors(a,:),'linewidth',1.5)

    end
    xlabel('Peak latency (ms)');
    ylabel('Cumulative count');
    formatForLee(gcf);
    set(gca,'fontsize',14)
    l=legend('5','10','15','20','25','30','40','50','100');
    set(l,'box','off','fontsize',14);
    
%% get number of spikes (probability) as a function of polarity and pulse width
% that for each neuron (on one plot probably)

    optsSpikesPlot.PRE_WINDOW = [-100,-5]/1000; % in s
    optsSpikesPlot.POST_WINDOW = [1,6]/1000; % in s
    optsSpikesPlot.PW1 = [200];
    optsSpikesPlot.PW2 = [200];
    optsSpikesPlot.POL = []; % 0 is cathodic first
    optsSpikesPlot.AMP = 50;
    optsSpikesPlot.PLOT_PULSE_WIDTH = 1;
    spikesStruct = {};
    spike_count = [];
    
    for unit = 1:numel(arrayData)
        if(unit == 1)
            optsSpikesPlot.MAKE_FIGURE = 1;
        else
            optsSpikesPlot.MAKE_FIGURE = 0;
        end
        [spikesStruct{unit},figure_handles] = plotSpikesPerCondition(arrayData{unit},optsSpikesPlot);
        hold on
        
        formatForLee(gcf)
        xlabel('First phase pulse width (\mus)');
        ylabel('Spikes per stimulation');
        set(gca,'fontsize',14)
        ax = gca;
       
        if(sum(spikesStruct{unit}.keep_mask) == 2) % should be anodic and cathodic entry
            spikes_struct_idx = find(spikesStruct{unit}.keep_mask);
            spike_count(end+1,:) = [sum(spikesStruct{unit}.spikes_post_stim{spikes_struct_idx(1)})/numel(spikesStruct{unit}.spikes_post_stim{spikes_struct_idx(1)}),...
                sum(spikesStruct{unit}.spikes_post_stim{spikes_struct_idx(2)})/numel(spikesStruct{unit}.spikes_post_stim{spikes_struct_idx(2)})];
        else
            spike_count(end+1,:) = [nan,nan];
        end
            
    end
    
    spike_count(:,3) = diff(spike_count')';
    
%% get number of spikes across the whole population at different latencies using a gaussian filter
% take advantage of filteredPSTH in plotLatencyActivation

    optsLatencyPlot.PW1 = 200;
    optsLatencyPlot.PW2 = 200;
    optsLatencyPlot.POL = 0; % 0 is cathodic first
    optsLatencyPlot.BIN_SIZE = 0.01; % ms
    optsLatencyPlot.KERNEL_LENGTH = 20; % in bins
    optsLatencyPlot.PEAK_WINDOW = [0,20]; % in ms
    optsLatencyPlot.BASELINE_WINDOW = [-100,-5]; % in ms
    
    optsLatencyPlot.ADJUST_LATENCY_TIME = 1;
    
    colors = inferno(numel(amps_plot));
    latencyStruct = {};
    
    amps_plot = [5,10,15,20,25,30,40,50,100];
    
    filtered_PSTH = [];
    num_units = zeros(numel(amps_plot),1);
    f_latency = figure(); hold on;
    
    for unit = 1:numel(arrayData)
        if(unit == 1)
            optsLatencyPlot.MAKE_FIGURE = 1;
        else
            optsLatencyPlot.MAKE_FIGURE = 0;
        end
        [latencyStruct{unit},figure_handles] = plotLatencyActivation(arrayData{unit},optsLatencyPlot);
        
        if(unit == 1) % initialize filtered_PSTH
            filtered_PSTH = zeros(numel(amps_plot),size(latencyStruct{unit}.filtered_PSTH,2));
        end
        
        for a = 1:numel(latencyStruct{unit}.amp)
            filtered_PSTH_idx = find(latencyStruct{unit}.amp(a) == amps_plot);

            filtered_PSTH(filtered_PSTH_idx,:) = filtered_PSTH(filtered_PSTH_idx,:) + latencyStruct{unit}.filtered_PSTH(a,:);
            num_units(filtered_PSTH_idx) = num_units(filtered_PSTH_idx) + 1;
        end
    end
    for a = 1:size(filtered_PSTH,1)
        plot(filtered_PSTH(a,:)/num_units(a),'color',colors(a,:),'linewidth',1.5);
        hold on
    end
    l=legend('5','10','15','20','25','30','40','50','100');
    
%% gaussian filter spike result to look at peaks specifically, look at independence of those peaks as well

    optsLatencyPlot.PW1 = 200;
    optsLatencyPlot.PW2 = 200;
    optsLatencyPlot.POL = 0; % 0 is cathodic first
    optsLatencyPlot.BIN_SIZE = 0.01; % ms
    optsLatencyPlot.KERNEL_LENGTH = 10; % in bins
    optsLatencyPlot.PEAK_WINDOW = [0,20]; % in ms
    optsLatencyPlot.BASELINE_WINDOW = [-100,-5]; % in ms
    
    optsLatencyPlot.ADJUST_LATENCY_TIME = 1;
    
    optsIndependence.PEAK_WIDTH = 0.5; % ms
    
    latencyStruct = {};
    peakIndependenceStruct = {};
    f_latency = figure(); hold on;
    f_peak_independence = figure();
    
    for unit = 1:numel(arrayData)
        if(unit == 1)
            optsLatencyPlot.MAKE_FIGURE = 1;
        else
            optsLatencyPlot.MAKE_FIGURE = 0;
        end
        [latencyStruct{unit},figure_handles] = plotLatencyActivation(arrayData{unit},optsLatencyPlot);
        [peakIndependenceStruct{unit}] = getPeakIndependence(arrayData{unit},latencyStruct{unit},optsIndependence);
        
        for a = 1:numel(latencyStruct{unit}.keep_mask)
            if(latencyStruct{unit}.keep_mask(a) == 1)
                plot(f_latency.Children,latencyStruct{unit}.amp(a)+zeros(size(latencyStruct{unit}.latencies{a})),...
                    latencyStruct{unit}.latencies{a},'.','markersize',20,'color','k')
                hold on
            end
        end
                
        formatForLee(f_latency)
        xlabel('Amplitude (\muA)');
        ylabel('Peak latencies (ms)');
        set(f_latency.Children,'fontsize',14)
        
        % plot independence calculation
        for a = 1:numel(latencyStruct{unit}.keep_mask)
            if(latencyStruct{unit}.keep_mask(a) == 1 && ~isempty(peakIndependenceStruct{unit}.pair_data{a}))
                plot(f_peak_independence.Children,peakIndependenceStruct{unit}.pair_data{a}(:,1),...
                    peakIndependenceStruct{unit}.pair_data{a}(:,2),'.','markersize',20,'color','k')
                hold on
            end
        end
       
        formatForLee(f_peak_independence)
        xlim([0,1])
        ylim([0,1])
        set(f_latency.Children,'fontsize',14)
        
    end


    
% make peak latency distribution plot -- must have latencyStruct (previous section)
    amps_plot = [5,10,15,20,25,30,40,50,100];
    pw1 = 200;
    pol = 0;
    latency_per_amp = {};
    bin_counts_per_amp = [];
    bin_edges = [0:2:20];
    colors = inferno(numel(amps_plot) + 3);
    f_latency=figure(); hold on;
    
    
    for a = 1:numel(amps_plot)
        latency_per_amp{a} = [];
        for unit = 1:numel(latencyStruct)
            amp_idx = find(latencyStruct{unit}.amp == amps_plot(a) & latencyStruct{unit}.pw1 == pw1 & latencyStruct{unit}.pol == pol);
            if(~isempty(amp_idx))
                latency_per_amp{a} = [latency_per_amp{end}; latencyStruct{unit}.latencies{amp_idx}'];
            end
        end
        [bin_counts_per_amp(:,a)] = histcounts(latency_per_amp{end},bin_edges);
        plot(bin_edges(1:end-1)+mode(diff(bin_edges)/2),bin_counts_per_amp(:,a),'color',colors(a,:),'linewidth',1.5)
        
    end
    
    xlabel('Peak latency after stimulation offset (ms)');
    ylabel('Count');
    formatForLee(gcf)
    set(gca,'fontsize',14)
    
% make cdf plot for latencies (dependent on previous section)
    amps_plot = [5,10,15,20,25,30,40,50,100];
    x_data = {};
    y_data = {};
    colors = inferno(numel(amps_plot) + 3);
    figure(); hold on
    for a = 1:numel(amps_plot)
        
        sorted_latencies = sort(latency_per_amp{a});
        
        x_data{a} = sort([0;latency_per_amp{a};latency_per_amp{a};max(bin_edges)]);
        y_data_list = (0:numel(latency_per_amp{a}));
        y_data{a} = reshape(repmat(y_data_list,2,1),[1,numel(y_data_list)*2]);
        
        plot(x_data{a},y_data{a},'color',colors(a,:),'linewidth',1.5)

    end
    xlabel('Peak latency (ms)');
    ylabel('Cumulative count');
    formatForLee(gcf);
    set(gca,'fontsize',14)
    
    
    
%% metric for inhibition
    optsInhibPlot = [];
    optsInhibPlot.PRE_WINDOW = [-50,-5];
    optsInhibPlot.POST_WINDOW = [0,200];
    optsInhibPlot.MAX_TIME_START = 40; % ms
    optsInhibPlot.BIN_SIZE = 1;
    optsInhibPlot.KERNEL_LENGTH = 5;
    optsInhibPlot.BLANK_TIME = [0,5]; % ms
    
    optsInhibPlot.PW1 = 200;
    optsInhibPlot.PW2 = 200;
    optsInhibPlot.POL = 0; % 0 is cathodic first
    optsInhibPlot.NUM_CONSECUTIVE_BINS = 10;
    
    inhibStruct = {};
    for unit = 1:numel(arrayData)
        if(unit == 1)
            optsSpikesPlot.MAKE_FIGURE = 1;
        else
            optsSpikesPlot.MAKE_FIGURE = 0;
        end
        [inhibStruct{unit},figure_handles] = plotInhibitionDuration(arrayData{unit},optsInhibPlot);

    end
    
    
    amp = [];
    unit = [];
    inhib_dur = [];
    for u = 1:numel(inhibStruct)
        amp = [amp,inhibStruct{u}.amp];
        inhib_dur = [inhib_dur,inhibStruct{u}.inhib_dur'];
        unit = [unit,u*ones(size(inhibStruct{u}.inhib_dur))'];
        
    end
    % remove nan's
    keep_mask = ~isnan(amp) & ~isnan(inhib_dur) & ~isnan(unit);
    amp = amp(keep_mask);
    inhib_dur = inhib_dur(keep_mask);
    unit = unit(keep_mask);
    
    tbl = table(amp',inhib_dur',unit','VariableNames',{'amp','inhib_dur','unit'});
    mdl = fitlm(tbl,'inhib_dur~amp+unit')
    
%% plot each unit as a faded out line, then plot mean as a dark line
    
    figure();
    faded_color = [0,0,0,0.5];
    amps_plot = [5,10,15,20,25,30,40,50,100];
    num_units_per_amp = zeros(size(amps_plot));
    inhib_dur_total = zeros(size(amps_plot));
    num_units_inhib = zeros(size(amps_plot));
    
    for unit = 1:numel(arrayData)
       for a = 1:numel(inhibStruct{unit}.amp)
            if(inhibStruct{unit}.keep_mask(a) == 1 && ~isempty(find(amps_plot == inhibStruct{unit}.amp(a))) && inhibStruct{unit}.is_inhib(a))
                amp_idx = find(amps_plot == inhibStruct{unit}.amp(a));
                inhib_dur_total(amp_idx) = inhib_dur_total(amp_idx) + inhibStruct{unit}.inhib_dur(a);
                
                
                num_units_inhib(amp_idx) = num_units_inhib(amp_idx) + 1;
            end
            
            if((inhibStruct{unit}.keep_mask(a) == 1 || (inhibStruct{unit}.amp(a) ~=50 && inhibStruct{unit}.amp(a) ~= 20)) && ~isempty(find(amps_plot == inhibStruct{unit}.amp(a))))
                amp_idx = find(amps_plot == inhibStruct{unit}.amp(a));
                num_units_per_amp(amp_idx) = num_units_per_amp(amp_idx) + 1;
            end
            
        end
        
        plot(inhibStruct{unit}.amp(inhibStruct{unit}.keep_mask==1),inhibStruct{unit}.inhib_dur(inhibStruct{unit}.keep_mask==1),...
            'marker','none','markersize',16,'color',faded_color,'linewidth',0.5)
        
        hold on

    end
    plot(amps_plot,inhib_dur_total./num_units_inhib,'k-','linewidth',2);
    formatForLee(gcf)
    xlabel('Amplitude (\muA)');
    ylabel('Inhibition duration (ms)');
    set(gca,'fontsize',14)
       
    
    figure();
    plot(amps_plot,num_units_inhib./num_units_per_amp,'marker','.','markersize',16,'linewidth',1.5)
    
    formatForLee(gcf)
    xlabel('Amplitude (\muA)');
    ylabel('Percent of cells with inhibitory response');
    set(gca,'fontsize',14)
       

%% 

    unit_idx = 25;
    waves = [9];
    figure()
    for w = waves
        plot(inhibStruct{unit}.filtered_PSTH(w,:))
        hold on
    end
    
    
    
    
    















