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

    for arrIdx = 1%:numel(arrayData)
%     arrIdx = 7;

        % plot raster, and PSTH for the given unit above
        
    %     optsPlotFunc.BIN_SIZE = optsExtract.BIN_SIZE;
        optsPlotFunc.BIN_SIZE = mode(diff(arrayData{arrIdx}.binEdges{1,1}));
        optsPlotFunc.FIGURE_SAVE = 0;
        optsPlotFunc.FIGURE_DIR = [];
        optsPlotFunc.FIGURE_PREFIX = [arrayData{arrIdx}.monkey,'_short_'];

        optsPlotFunc.PRE_TIME = 30/1000;
        optsPlotFunc.POST_TIME = 100/1000;
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

    end


%% get number of spikes (probability) as a function of amplitude and plot
% that for each neuron (on one plot probably)

    optsSpikesPlot.PRE_WINDOW = [-81,-5]/1000; % in s
    optsSpikesPlot.POST_WINDOW = [0,5]/1000; % in s, for amp plot
%     optsSpikesPlot.POST_WINDOW = [1,5]/1000; % in s, for latency plot
    optsSpikesPlot.MAXIMIZE_RESPONSE_PROB = 0; % find maximum for each condition instead of using the window
%     optsSpikesPlot.POST_WINDOW = [0.75,3; % 5uA
%                                     0.75,3; % 10uA
%                                     0.75,3; % 15uA
%                                     0.75,5; % 20uA
%                                     0.75,5; % 25uA
%                                     0.75,7; % 30uA
%                                     1,7.5; % 40uA
%                                     1.25,7.5; % 50uA
%                                     1.25,7.5; % 60uA
%                                     1.5,8.25; % 80uA
%                                     1.5,9]/1000; % 100uA
%     optsSpikesPlot.POST_WINDOW(:,1) = 1.5/1000;
    optsSpikesPlot.AMP_LIST = [5;10;15;20;25;30;40;50;60;80;100];
    
    optsSpikesPlot.PW1 = 200;
    optsSpikesPlot.PW2 = 200;
    optsSpikesPlot.POL = 0; % 0 is cathodic first
    optsSpikesPlot.AMP = [];
    optsSpikesPlot.PLOT_PULSE_WIDTH = 0;
    
    optsSpikesPlot.ADJUST_SPIKE_TIMES = 1;
    optsSpikesPlot.WAVEFORM_LENGTH = 0.453/1000; % s
    
    spikesStruct = {}; amp = []; spikes_evoked = []; unit_idx = [];
    max_rgb = 0.6;
    min_rgb = 0;
    
    for unit = 1:numel(arrayData)
        if(unit == 1)
            optsSpikesPlot.MAKE_FIGURE = 1;
        else
            optsSpikesPlot.MAKE_FIGURE = 0;
        end
        optsSpikesPlot.COLOR = ones(1,3)*(unit/numel(arrayData))*(max_rgb-min_rgb) + min_rgb; 
        if(unit == 15)
            optsSpikesPlot.COLOR = getColorFromList(1,1);
        end
        [spikesStruct{unit},figure_handles] = plotSpikesPerCondition(arrayData{unit},optsSpikesPlot);
        hold on
        
        formatForLee(gcf)
        xlabel('Amplitude (\muA)');
        ylabel('Spikes per stimulation');
        set(gca,'fontsize',14)
        
        amp = [amp,spikesStruct{unit}.amp(spikesStruct{unit}.keep_mask==1)];
        spikes_evoked = [spikes_evoked, spikesStruct{unit}.prob_spike(spikesStruct{unit}.keep_mask==1)];
        unit_idx = [unit_idx, unit + zeros(1,sum(sum(spikesStruct{unit}.keep_mask==1)))];
    end  
%     
    amps_plot = unique(amp);
    mean_spikes_evoked = zeros(1,numel(amps_plot));
    
    for a = 1:numel(amps_plot)
        amp_mask = amp == amps_plot(a);
        mean_spikes_evoked(a) = mean(spikes_evoked(amp_mask));
        
    end
    
    tbl = table(amp',spikes_evoked',unit_idx','VariableNames',{'amp','spikes_evoked','unit'});
    tbl.unit = categorical(tbl.unit);
    mdl = fitlm(tbl,'spikes_evoked~amp+unit')
    
    x_data = unique(amp);
    y_data = mean(mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(3:end)) + mdl.Coefficients.Estimate(2).*x_data;
    plot(x_data,y_data,'color',getColorFromList(1,0),'linewidth',3,'linestyle','--');
%     figure();
%     plot(mdl)
    
%% compare response at two amplitudes
    amps = [25,50];
    data = [];
    for unit = 1:numel(spikesStruct)
        amp_1_idx = find(spikesStruct{unit}.amp == amps(1) & spikesStruct{unit}.keep_mask == 1);
        amp_2_idx = find(spikesStruct{unit}.amp == amps(2) & spikesStruct{unit}.keep_mask == 1);
        
        if(~isempty(amp_1_idx) && ~isempty(amp_2_idx))
            data = [data; spikesStruct{unit}.prob_spike(amp_1_idx), spikesStruct{unit}.prob_spike(amp_2_idx)];
        end
    end
    
    figure();
    subplot(1,2,1)
    plot(data','k')
    subplot(1,2,2)
    histogram(diff(data'),-2:0.2:2);
    
%% find threshold (traditional 50% definition)
    threshold = nan(numel(spikesStruct),1);
    counter = 0;
    for i_unit = 1:numel(spikesStruct)
        percent_respond = cellfun(@mean,spikesStruct{i_unit}.num_spikes_post_stim);
        if(~isempty(find(percent_respond >= 0.5)))
            counter = counter + 1;
            threshold(i_unit) = spikesStruct{i_unit}.amp(find(percent_respond >= 0.5,1,'first'));
        end
    end

%% make spike time distribution plot, and % excitatory response at different amplitudes
    amps_plot = [5,10,15,20,25,30,40,50,100];
    spike_times_amp = cell(numel(amps_plot),1); % preallocate space  
    num_responsive = zeros(1,numel(amps_plot));
    num_total = zeros(1,numel(amps_plot));
    num_stims = zeros(1,numel(amps_plot));
    
    bin_edges = [0:1:20];
    bin_counts = cell(numel(amps_plot),1);
    colors = inferno(10);
    f_latency=figure(); hold on;
    simulation_gain = ones(numel(amps_plot),size(bin_edges,2)-1);
    if(exist('spikes_found'))
        [~,~,bin_edge_idx] = histcounts(spikes_found_t,bin_edges);
        unique_bin_edge_idx = unique(bin_edge_idx);
        unique_bin_edge_idx(unique_bin_edge_idx == 0) = [];
        for i_amp = 1:numel(amps_plot)
            for i_unique = unique_bin_edge_idx
%                 simulation_gain(i_amp,i_unique) = 1/mean(spikes_found(i_amp,bin_edge_idx == i_unique)./spikes_total(i_amp,bin_edge_idx == i_unique),'omitnan');

            end
        end
    end
    
    simulation_gain(isnan(simulation_gain)) = 1;
    color_counter= 1;
    baseline_fr_all = cell(numel(amps_plot),1);
    for a = 1:numel(amps_plot)
        bin_counts{a} = zeros(1,numel(bin_edges)-1);
        for unit = 1:numel(spikesStruct)
            amp_idx = find(spikesStruct{unit}.amp == amps_plot(a) & spikesStruct{unit}.keep_mask == 1);
            unit_param_idx = find([array_data_all{unit}.STIM_PARAMETERS.amp1] == amps_plot(a));

            if(~isempty(amp_idx) && ~isempty(unit_param_idx))
                num_total(a) = num_total(a) + 1;
                if(~isnan(threshold(unit)) && threshold(unit) < 100)
%                 if(any(spikesStruct{unit}.is_excitatory_p < 0.05) || 1==1) % use all neurons for now....
                    bin_counts{a} = bin_counts{a} + histcounts(spikesStruct{unit}.spike_times_post_stim{amp_idx}*1000,bin_edges)/array_data_all{unit}.numStims(unit_param_idx(1));
                    spike_times_amp{a} = [spike_times_amp{a},spikesStruct{unit}.spike_times_post_stim{amp_idx}*1000];
                    num_responsive(a) = num_responsive(a) + 1;% (spikesStruct{unit}.is_excitatory_p(amp_idx) < 0.05);
                    num_stims(a) = num_stims(a) + arrayData{unit}.numStims(amp_idx);   
                        baseline_fr_all{a}(end+1,1) = spikesStruct{unit}.mean_baseline_fr;
                end
            end
        end
        plot(bin_edges(1:end-1)+mode(diff(bin_edges)/2),bin_counts{a}.*simulation_gain(a,:)/num_responsive(a),'color',colors(color_counter,:),'linewidth',1.5)
%         plot(simulation_gain(a,:),'color',colors(a,:),'linewidth',1.5)
        hold on
        color_counter = color_counter + 1;
    end
    
    xlabel('Time after stimulation offset (ms)');
    ylabel({'Number of spikes per stimulation','per neuron'});
    formatForLee(gcf)
    set(gca,'fontsize',14)
%     l=legend('5','10','15','20','25','30','40','50','100');
    l=legend('5\muA','10\muA','15\muA','25\muA','40\muA','100\muA')
    set(l,'box','off','fontsize',14);
    
%     f_percent_responsive = figure();
%     plot(amps_plot,num_responsive./num_total)
    
%% make cdf plot for latencies (dependent on previous section)
    
    x_data = {};
    y_data = {};
    colors = inferno(10);
    figure(); hold on
    color_counter = 1;
    linewidth_list = linspace(0.75,3,10);
    for a = 1:numel(amps_plot)
                
        x_data{a} = sort([0;spike_times_amp{a}';spike_times_amp{a}';max(bin_edges)]);
        y_data_list = (0:numel(spike_times_amp{a}))/num_stims(a);
        y_data{a} = reshape(repmat(y_data_list,2,1),[1,numel(y_data_list)*2]);
        
        plot(x_data{a},y_data{a},'color',colors(color_counter,:),'linewidth',linewidth_list(color_counter))
        color_counter = color_counter + 1;
    end
    xlabel('Time after stim offset (ms)');
    ylabel('Cumulative spikes per stimulation');
    formatForLee(gcf);
    set(gca,'fontsize',14)
%     l=legend('5','10','15','20','25','30','40','50','100');
    l=legend('5\muA','10\muA','15\muA','25\muA','40\muA','100\muA')
    set(l,'box','off','fontsize',14);
    xlim([0,10])
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
    optsLatencyPlot.KERNEL_LENGTH = 50; % in bins
    optsLatencyPlot.PEAK_WINDOW = [0,20]; % in ms
    optsLatencyPlot.BASELINE_WINDOW = [-80,-5]; % in ms
    
    optsLatencyPlot.ADJUST_LATENCY_TIME = 1;
    
    colors = inferno(numel(amps_plot)+2);
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
        plot(latencyStruct{1}.x_data-0.453, filtered_PSTH(a,:)/num_units(a),'color',colors(a,:),'linewidth',1.5);
        hold on
    end
    l=legend('5','10','15','20','25','30','40','50','100');
    xlim([0,20])

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
    optsInhibPlot.PRE_WINDOW = [-80,-5];
    optsInhibPlot.POST_WINDOW = [0,220];
    optsInhibPlot.MAX_TIME_START = 40; % ms
    optsInhibPlot.BIN_SIZE = 1;
    optsInhibPlot.KERNEL_LENGTH = 10;
    optsInhibPlot.BLANK_TIME = [0,5]; % ms
    
    optsInhibPlot.PW1 = 200;
    optsInhibPlot.PW2 = 200;
    optsInhibPlot.POL = 0; % 0 is cathodic first
    optsInhibPlot.NUM_CONSECUTIVE_BINS = 10;
    optsInhibPlot.AMPS_PLOT = [5,10,15,20,25,30,40,50,100];
    
    inhibStruct = {};
    for unit = 1:104%numel(arrayData)
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
%     % remove nan's
    keep_mask = ~isnan(amp) & ~isnan(inhib_dur) & ~isnan(unit);
    amp = amp(keep_mask);
    inhib_dur = inhib_dur(keep_mask);
    unit_idx = unit(keep_mask);
%     
% %     [inhib_dur_adj,lambda] = boxcox(inhib_dur');
%     
    inhib_dur_all = inhib_dur;
    % plot inhib dur vs. baseline fr for a single amp
%     amp_plot = 100;
%     inhib_dur = [];
%     baseline_fr = [];
%     for u = 1:numel(inhibStruct)
%         amp_idx = find(inhibStruct{u}.amp == amp_plot,1,'first');
%         if(~isempty(amp_idx) && inhibStruct{u}.is_inhib(amp_idx)==1)
%             baseline_fr(end+1,1) = inhibStruct{u}.baseline_fr;
%             inhib_dur(end+1,1) = inhibStruct{u}.inhib_dur(amp_idx);
%         end
%     end
    
%     figure(); plot(baseline_fr,inhib_dur,'.','markersize',20)
% plot each unit as a faded out line, then plot mean as a dark line
    max_rgb = 0.6;
    min_rgb = 0;
    num_colors = numel(arrayData);
    f=figure(); f.Name = 'Duncan_Han_inhib_dur';
    
    amps_plot = [5,10,15,20,25,30,40,50,100];
    num_units_per_amp = zeros(2,numel(amps_plot));
    inhib_dur_total = zeros(size(amps_plot));
    num_units_inhib = zeros(2,numel(amps_plot));
    monkey = ones(numel(arrayData),1);
    
    for unit = 1:numel(arrayData)
        if(strcmpi(arrayData{unit}.monkey,'Duncan'))
            monkey(unit) = 2;
        end
        for a = 1:numel(inhibStruct{unit}.amp)
            if(inhibStruct{unit}.keep_mask(a) == 1 && ~isempty(find(amps_plot == inhibStruct{unit}.amp(a))) && inhibStruct{unit}.is_inhib(a))
                amp_idx = find(amps_plot == inhibStruct{unit}.amp(a));
                inhib_dur_total(amp_idx) = inhib_dur_total(amp_idx) + inhibStruct{unit}.inhib_dur(a);

                num_units_inhib(monkey(unit),amp_idx) = num_units_inhib(monkey(unit),amp_idx) + 1;
            end
            
            if((inhibStruct{unit}.keep_mask(a) == 1 || (inhibStruct{unit}.amp(a) ~=50 && inhibStruct{unit}.amp(a) ~= 20)) && ~isempty(find(amps_plot == inhibStruct{unit}.amp(a))))
                amp_idx = find(amps_plot == inhibStruct{unit}.amp(a));
                num_units_per_amp(monkey(unit),amp_idx) = num_units_per_amp(monkey(unit),amp_idx) + 1;
            end
            
        end
        faded_color = ones(1,3).*(unit/num_colors)*(max_rgb-min_rgb) + min_rgb;
        if(unit == 15)
            faded_color = getColorFromList(1,1);
        end
        plot(inhibStruct{unit}.amp(inhibStruct{unit}.keep_mask==1),inhibStruct{unit}.inhib_dur(inhibStruct{unit}.keep_mask==1),...
            'marker','.','markersize',12,'color',faded_color,'linewidth',0.75)
        hold on

    end
    
    tbl = table(amp',sqrt(inhib_dur_all'),unit_idx','VariableNames',{'amp','inhib_dur','unit'});
    tbl.unit = categorical(tbl.unit);
    mdl = fitlm(tbl,'inhib_dur~amp+unit')
    x_data = unique(amp);
    y_data = mean(mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(3:end)) + mdl.Coefficients.Estimate(2).*x_data;
    plot(x_data,y_data.^2,'color',getColorFromList(1,0),'linewidth',3,'linestyle','--');
    
    formatForLee(gcf)
    xlabel('Amplitude (\muA)');
    ylabel('Inhibition duration (ms)');
    set(gca,'fontsize',14)
       
    
    f=figure(); hold on;
    f.Name = 'Duncan_Han_inhibitory_percent_cell';
%     for i = 1:2
        plot(amps_plot,sum(num_units_inhib,1)./sum(num_units_per_amp,1),'marker','.','markersize',16,'linewidth',1.5,'linestyle','-','color',getColorFromList(1,1))
%         fits = fit(amps_plot',(num_units_inhib(i,:)./num_units_per_amp(i,:))','cubicinterp')
%         plot(amps_plot,feval(fits,amps_plot),'linewidth',1.5,'color',getColorFromList(1,1+i))
%     end
    formatForLee(gcf)
    xlabel('Amplitude (\muA)');
    ylabel('Percent of cells with inhibitory response');
    set(gca,'fontsize',14)
%     l=legend('Monkey H','Monkey D');
%     set(l,'box','off')
    ylim([0,1])
    xlim([0,105])
%% plot PSTH for each condition for low and high speed cases
    mean_speed_cutoff = [0.2,4]; % plot speeds in range
    input_data.suffix = 'speed0-2';
    
    array_data_trim = arrayData;
    for u = 1:numel(arrayData)
        for cond = 1:numel(arrayData{u}.binCounts)
            arrayData{u}.kin{cond}.speed = sqrt(arrayData{u}.kin{cond}.vx.^2 + arrayData{u}.kin{cond}.vy.^2);
            arrayData{u}.kin{cond}.mean_speed = mean(arrayData{u}.kin{cond}.speed,2);
            trial_speeds = arrayData{u}.kin{cond}.mean_speed(arrayData{u}.stimData{cond});
            keep_mask = trial_speeds > mean_speed_cutoff(1) & trial_speeds < mean_speed_cutoff(2);
            trials_keep = unique(arrayData{u}.stimData{cond}(keep_mask));
            
            
            spike_trial_times = arrayData{u}.spikeTrialTimes{cond}(keep_mask==1);
            
            array_data_trim{u}.num_stims(cond) = sum(arrayData{u}.kin{cond}.mean_speed > mean_speed_cutoff(1) & arrayData{u}.kin{cond}.mean_speed < mean_speed_cutoff(2));
            array_data_trim{u}.binCounts{cond} = histcounts(spike_trial_times*1000,array_data_trim{u}.binEdges{cond})/array_data_trim{u}.num_stims(cond);
            array_data_trim{u}.speed{cond} = arrayData{u}.kin{cond}.speed(trials_keep,:);
            array_data_trim{u}.mean_speed{cond} = arrayData{u}.kin{cond}.mean_speed(trials_keep,:);
            
            array_data_trim{u}.binMaxYLim = 1;
        end
    end
    
    input_data.unit_idx = 1;
    
    optsPlotFunc.BIN_SIZE = mode(diff(array_data_trim{arrIdx}.binEdges{1,1}));
    optsPlotFunc.FIGURE_SAVE = 0;
    optsPlotFunc.FIGURE_DIR = [];
    optsPlotFunc.FIGURE_PREFIX = [array_data_trim{arrIdx}.monkey,'_long_'];

    optsPlotFunc.PRE_TIME = 15/1000;
    optsPlotFunc.POST_TIME = 20/1000;
    optsPlotFunc.SORT_DATA = 'postStimuliTime';

    optsPlotFunc.MARKER_STYLE  = '.';

    optsPlotFunc.PLOT_AFTER_STIMULATION_END = 1;
    optsPlotFunc.STIMULATION_LENGTH = [];

    optsPlotFunc.PLOT_ALL_ONE_FIGURE = 0;
    optsPlotFunc.PLOT_LINE = 1;
    optsPlotFunc.PLOT_TITLE = 0;    
    optsPlotFunc.PLOT_ALL_WAVES_ONE_FIGURE = 0;
        
    %   
    PSTHPlots = plotPSTHStim(array_data_trim{arrIdx},array_data_trim{arrIdx}.NN,optsPlotFunc);    

    


%% plot response vs baseline firing rate for each trial

    num_spikes_per_trial = cell(numel(arrayData{1}.numStims),numel(arrayData));
    num_spikes_baseline = cell(numel(arrayData{1}.numStims),numel(arrayData));
    mean_speed_per_trial = cell(numel(arrayData{1}.numStims),numel(arrayData));
    spike_window = [0,10]/1000; % in s
    baseline_window = [-100,-5]/1000; % in s
    
    figure();
    plot([0,10],[0,10],'k--');
    hold on
    for unit = 1:numel(arrayData)
        for cond = 1:numel(arrayData{unit}.numStims)
            arrayData{u}.kin{cond}.speed = sqrt(arrayData{u}.kin{cond}.vx.^2 + arrayData{u}.kin{cond}.vy.^2);
            arrayData{u}.kin{cond}.mean_speed = mean(arrayData{u}.kin{cond}.speed,2);
            if(arrayData{unit}.numStims(cond) > 0)
                for tr = 1:arrayData{unit}.numStims(cond)
                    % get spikes during the train
                    spike_mask = arrayData{unit}.spikeTrialTimes{cond} > spike_window(1) & arrayData{unit}.spikeTrialTimes{cond} < spike_window(2) & ...
                        arrayData{unit}.stimData{cond} == tr;
                    num_spikes_per_trial{cond,unit}(tr,1) = sum(spike_mask);
                    mean_speed_per_trial{cond,unit}(tr,1) = arrayData{unit}.kin{cond}.mean_speed(tr);

                    % get baseline spike rate
                    spike_mask = arrayData{unit}.spikeTrialTimes{cond} > baseline_window(1) & arrayData{unit}.spikeTrialTimes{cond} < baseline_window(2) & ...
                        arrayData{unit}.stimData{cond} == tr;
                    num_spikes_baseline{cond,unit}(tr,1) = sum(spike_mask);  

                end
                
                % plot mean fr baseline for trials with >= 1 spike vs. that for stims with 0 spikes
                spike_trial_mask = num_spikes_per_trial{cond,unit} > 0;
                
                plot(mean(num_spikes_baseline{cond,unit}(spike_trial_mask==0)),mean(num_spikes_baseline{cond,unit}(spike_trial_mask==1)),'.','markersize',16)
            end
        end
    end

    




%% proportion of cells with only inhib response, only excite response, or both
%     is_excite_inhib = zeros(numel(spikesStruct),9,2); % 9 amps, is_excite, is_inhib
    is_excite_inhib = [];
    for u = 1:numel(spikesStruct)
        if(numel(inhibStruct{u}.amp) > 9)
            for a_idx = 1:size(spikesStruct{u}.is_excitatory_p,1)
                is_excite_inhib(end+1,:,1) = spikesStruct{u}.is_excitatory_p(a_idx,:) < 0.05;
                is_excite_inhib(end,:,2) = inhibStruct{u}.is_inhib(a_idx,:);
            end
        end
    end
    total_cells = size(is_excite_inhib,1);
    excite_only = sum(is_excite_inhib(:,:,1) == 1 & is_excite_inhib(:,:,2) == 0,1)./total_cells;
    inhib_only = sum(is_excite_inhib(:,:,1) == 0 & is_excite_inhib(:,:,2) == 1,1)./total_cells;
    excite_and_inhib = sum(is_excite_inhib(:,:,1) == 1 & is_excite_inhib(:,:,2) == 1,1)./total_cells;
    neither = sum(is_excite_inhib(:,:,1) == 0 & is_excite_inhib(:,:,2) == 0,1)./total_cells;
    
    x = inhibStruct{1}.amp(1:4:end);
    
    figure(); hold on
    
    plot(x(1:4),excite_only','color',getColorFromList(1,0),'linewidth',2)
    plot(x(1:4),inhib_only,'color',getColorFromList(1,1),'linewidth',2)
    plot(x(1:4),excite_and_inhib,'color',getColorFromList(1,2),'linewidth',2)
    plot(x(1:4),neither,'color',getColorFromList(1,3),'linewidth',2)
    l=legend('excite only','inhib only','both','neither');
    set(l,'box','off','location','best');
    set(gca,'fontsize',14);
    
    xlabel('Amplitude (\muA)')
    ylabel('Proportion')
    ylim([0,1])
    formatForLee(gcf);
%     x = repmat(x,3,1);
%     x = reshape(x,size(x,1)*size(x,2),1) + repmat([-0.3;0;0.3],size(is_excite_inhib,2),1);
    
    
    