%% set file names 

    input_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Han\Han_20190304_dblPulse\chan44stim\';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    % input_data.mapFileName = 'mapFileR:\limblab-archive\Retired Animal Logs\Monkeys\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
    % input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';


    input_data.IPI = [-1,10,20,200];
    input_data.num_conditions = numel(input_data.IPI);

    input_data.nom_freq = 2;
    input_data.window = [-50,450]; % 20ms before and 400ms after 1st of 2 pulses
    input_data.bin_size = 2; % in ms
    input_data.chan_rec = 44;

    folderpath = input_data.folderpath; % rest of code uses folderpath currently...may have switched this, not 100% certain

    input_data.task='taskCObump';
    input_data.ranBy='ranByJoseph'; 
    input_data.array1='arrayLeftS1'; 
    input_data.monkey='monkeyHan';
    input_data.labnum = 6;

    pwd=cd;
    cd(input_data.folderpath)
    fileList = dirSorted('*spikesExtracted.nev*');
    stimInfoFileList = dirSorted('*stimInfo*');


%% load in files, parse appropriately, then combine
    array_data = {};

    for fileNumber = 1:numel(fileList)
        disp(fileList(fileNumber).name)
        cd(input_data.folderpath)
        cds = commonDataStructure();
        load(stimInfoFileList(fileNumber).name);
        if(~iscell(stimInfo.chanSent))
            stimInfo.chanSent = mat2cell(stimInfo.chanSent',ones(numel(stimInfo.chanSent),1));
            warning('made chan sent a cell array');
        end
        cd(input_data.folderpath)
        cds.file2cds([input_data.folderpath fileList(fileNumber).name],input_data.task,input_data.ranBy,...
            input_data.monkey,input_data.labnum,input_data.array1,input_data.mapFileName); % DO NOT USE RECOVER PRE SYNC, currently this shifts the units and the analog signal differently

        if(fileNumber == 1)
            unit_idx = find([cds.units.chan] == input_data.chan_rec & [cds.units.ID] ~= 0 & [cds.units.ID] ~= 255);
            for u = 1:numel(unit_idx)
                array_data{u} = [];
                array_data{u}.spikeTrialTimes = cell(input_data.num_conditions,1);
                array_data{u}.trial_num = cell(input_data.num_conditions,1);
                array_data{u}.num_stims = zeros(input_data.num_conditions,1);
                array_data{u}.binCounts = cell(input_data.num_conditions,1);
                array_data{u}.binEdges = cell(input_data.num_conditions,1);
            end
        end
        % make stimOn a timestamp
        t = (0:1:stimInfo.stimOn(end))/30000;
        stimInfo.stimOn = t(stimInfo.stimOn);


        is_1st_pulse_mask = ones(numel(stimInfo.stimOn),1);

        for st = 1:numel(stimInfo.stimOn)
            if(is_1st_pulse_mask(st) == 1) % look for condition based on IPI to next stimulation
                if(st == numel(stimInfo.stimOn)) % last stim was a single pulse
                    condition = find(input_data.IPI < 0);
                else
                    IPI = 1000*(stimInfo.stimOn(st+1) - stimInfo.stimOn(st)); % in ms now, like input_data.IPI
                    condition = [];
                    [min_diff,condition] = min(abs(IPI - input_data.IPI));
                    if(min_diff > 50)
                        condition = find(input_data.IPI < 0);
                    end
                    if(input_data.IPI(condition) > 0) % double pulse case
                        is_1st_pulse_mask(st+1) = 0;
                    end
                end

                % we know which condition and have dealt with the double pulse
                % stuff, so now all we need to do is extract relevant spikes based
                % on input_data.window and store them into array_data

                for u = 1:numel(unit_idx)
                    spike_mask = cds.units(unit_idx(u)).spikes.ts > stimInfo.stimOn(st) + input_data.window(1)/1000 & cds.units(unit_idx(u)).spikes.ts < stimInfo.stimOn(st) + input_data.window(2)/1000;
                    spike_ts = cds.units(unit_idx(u)).spikes.ts(spike_mask == 1) - stimInfo.stimOn(st);
                    array_data{u}.spikeTrialTimes{condition} = [array_data{u}.spikeTrialTimes{condition}, spike_ts'];
                    array_data{u}.num_stims(condition) = array_data{u}.num_stims(condition) + 1;
                    array_data{u}.trial_num{condition} = [array_data{u}.trial_num{condition}, ones(1,numel(spike_ts))*array_data{u}.num_stims(condition)];
                end
            end

        end

        % get bin counts for each condition
        bin_edges = input_data.window(1):input_data.bin_size:input_data.window(2);

        for u = 1:numel(array_data)
            for cond = 1:numel(array_data{u}.binCounts)
                array_data{u}.binEdges{cond} = bin_edges;
                array_data{u}.binCounts{cond} = histcounts(array_data{u}.spikeTrialTimes{cond}*1000,bin_edges);
            end
        end
    end


%% rebin if necessary

    input_data.bin_size = 5; % ms
    bin_edges = input_data.window(1):input_data.bin_size:input_data.window(2);
    for u = 1:numel(array_data)
        for cond = 1:numel(array_data{u}.binCounts)
            array_data{u}.binEdges{cond} = bin_edges;
            array_data{u}.binCounts{cond} = histcounts(array_data{u}.spikeTrialTimes{cond}*1000,bin_edges);
        end
    end
    
     
%% plot raster 
    for arrIdx = 4%:numel(array_data)   

        optsPlotFunc.BIN_SIZE = mode(diff(array_data{arrIdx}.binEdges{1,1}));
        optsPlotFunc.FIGURE_SAVE = 0;
        optsPlotFunc.FIGURE_DIR = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\';
        optsPlotFunc.FIGURE_PREFIX = [array_data{arrIdx}.monkey,'_DblPulseTrains_'];

        optsPlotFunc.PRE_TIME = 75/1000;
        optsPlotFunc.POST_TIME = 300/1000;
        optsPlotFunc.SORT_DATA = '';

        optsPlotFunc.STIMULATION_LENGTH = 0.453;

        optsPlotFunc.MARKER_STYLE  = '.';

        optsPlotFunc.PLOT_AFTER_STIMULATION_END = 1;

        rasterPlots = plotRasterStim(array_data{arrIdx},arrIdx,optsPlotFunc);
    end
        
%% plot PSTH for each condition
    for unit_idx = 1:numel(array_data)
        input_data.unit_idx = unit_idx;
        input_data.chan_rec = array_data{unit_idx}.CHAN_LIST;
        plotPSTHArrayData(array_data,ifpnput_data);
    end
    
    
%% count spikes in window after each stimulation, plot those values
    baseline_window = [-100,-5]; % in ms
    window = [0,9]; % in ms
    post_stim_fr = nan(numel(array_data),input_data.num_conditions,2);
    stim_firing_rate = nan(numel(array_data),input_data.num_conditions,2);
    baseline_firing_rate = nan(numel(array_data),input_data.num_conditions,1);
    
    keep_mask = input_data.num_pulses == 2;
    single_pulse_idx = find(input_data.num_pulses == 1);
    monkey_idx = [];
    
    for u = 1:numel(array_data)
        for cond = 1:numel(array_data{u}.spikeTrialTimes)
            if(keep_mask(cond) == 1 || cond == single_pulse_idx)
                stim_offset = input_data.IPI(cond);
                if(input_data.IPI(cond) == 20)
                    stim_offset = 20.63;
                elseif(input_data.IPI(cond) == 10)
                    stim_offset = 10.6;
                elseif(input_data.IPI(cond) == 200)
                    stim_offset = 201.3;
                end
                window_adj = window + stim_offset;
                       
                baseline_firing_rate(u,cond) = sum(array_data{u}.spikeTrialTimes{cond}*1000 > baseline_window(1) & ...
                     array_data{u}.spikeTrialTimes{cond}*1000 < baseline_window(2))/ ...
                     array_data{u}.num_stims(cond)/(diff(baseline_window))*1000; % convert to Hz
                 
                stim_firing_rate(u,cond,1) = sum(array_data{u}.spikeTrialTimes{cond}*1000 > window(1) & ...
                     array_data{u}.spikeTrialTimes{cond}*1000 < window(2))/ ...
                     array_data{u}.num_stims(cond)/(diff(window))*1000; % convert to Hz, first pulse response
                
             
                stim_firing_rate(u,cond,2) = sum(array_data{u}.spikeTrialTimes{cond}*1000 > window_adj(1) & ...
                     array_data{u}.spikeTrialTimes{cond}*1000 < window_adj(2))/ ...
                     array_data{u}.num_stims(cond)/(diff(window_adj))*1000; % convert to Hz, second pulse response
                 
                post_stim_fr(u,cond,:) = stim_firing_rate(u,cond,:) - baseline_firing_rate(u,cond,1);
               
            end
        end
        
        monkey_idx(u) = strcmpi(array_data{u}.monkey,'Han'); % 1 = han, 0 = duncan
    end
    
    [sort_IPI,sort_idx] = sort(input_data.IPI);
    prob_spike_sort = post_stim_fr(:,sort_idx,:);
    keep_mask_sort = keep_mask(sort_idx);
  
%     figure();
%     plot(sort_IPI(keep_mask_sort == 1),prob_spike_sort(:,keep_mask_sort == 1,2) - prob_spike_sort(:,keep_mask_sort == 1,1))

% dot plot response to single pulse vs. second pulse at different latencies
    markers = {'.','s'};
    marker_size = [22,8];
    f=figure();
    f.Name = 'Duncan_Han_dblPulse_excitation_summary';
    hold on
    for m = 1:2 % duncan, then han
        keep_mask = monkey_idx == m-1; % matlab indexing causes issues...
        
        plot(post_stim_fr(keep_mask,single_pulse_idx,1),post_stim_fr(keep_mask,input_data.num_pulses == 2 & input_data.IPI == 200,2),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,0),'markerfacecolor',getColorFromList(1,0),'markersize',marker_size(m));
        
        plot(post_stim_fr(keep_mask,single_pulse_idx,1),post_stim_fr(keep_mask,input_data.num_pulses == 2 & input_data.IPI == 20,2),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,1),'markerfacecolor',getColorFromList(1,1),'markersize',marker_size(m));
        
        plot(post_stim_fr(keep_mask,single_pulse_idx,1),post_stim_fr(keep_mask,input_data.num_pulses == 2 & input_data.IPI == 10,2),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,2),'markerfacecolor',getColorFromList(1,2),'markersize',marker_size(m));
        
    end
    
    unity_line = plot([-20,150],[-20,150],'k--','linewidth',1.5);

    xlabel('Response to single pulse (Hz)');
    ylabel('Response to second pulse (Hz)');
    l=legend('200ms','20ms','10ms'); % actual IPI is [10.6,20.6,201.03]
    set(l,'box','off','fontsize',14,'location','northwest');
    uistack(unity_line,'bottom');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    xlim([-20,150]); ylim([-20,150])
%     ylim([0,1]);


%% look at duration of inhibition after second pulse compared to the single pulse case
    optsInhibPlot = [];
    optsInhibPlot.PRE_WINDOW = [-200,-5];
    post_window = [0,200]; % ms
    optsInhibPlot.MAX_TIME_START = 70; % ms
    optsInhibPlot.BIN_SIZE = 1;
    optsInhibPlot.KERNEL_LENGTH = 5;
    blank_time = [0,5]; % ms
    
    optsInhibPlot.PW1 = 200;
    optsInhibPlot.PW2 = 200;
    optsInhibPlot.POL = 0; % 0 is cathodic first
    optsInhibPlot.NUM_CONSECUTIVE_BINS = 10;
    
    inhibStruct = {};
    keep_mask = input_data.num_pulses <= 2;
    
    figure();
    inhib_dur = [];
    monkey_idx = [];
    for u = 1:numel(array_data)
        for cond = 1:numel(array_data{u}.spikeTrialTimes)
            if(input_data.num_pulses(cond) == 1) % single pulse
                optsInhibPlot.POST_WINDOW(cond,:) = post_window;
                optsInhibPlot.BLANK_TIME(cond,:) = blank_time;
            elseif(input_data.num_pulses(cond) == 2)% double pulse
                optsInhibPlot.POST_WINDOW(cond,:) = post_window + input_data.IPI(cond);
                optsInhibPlot.BLANK_TIME(cond,:) = blank_time + input_data.IPI(cond);
            else % train
                optsInhibPlot.POST_WINDOW(cond,:) = post_window + 200;
                optsInhibPlot.BLANK_TIME(cond,:) = blank_time + 200;

            end
        
        end
        
        [inhibStruct{u},figure_handles] = plotInhibitionDuration(array_data{u},optsInhibPlot);

        monkey_idx(u) = strcmpi(array_data{u}.monkey,'Han'); % 1 = han, 0 = duncan 
        inhib_dur(end+1,:) = inhibStruct{u}.inhib_dur;
        
        % plot inhib duration for each condition, have to sort IPI first
        [IPI_sort,sort_idx] = sort(input_data.IPI);
        inhib_dur_sort = inhibStruct{u}.inhib_dur(sort_idx);
        keep_mask_sort = keep_mask(sort_idx);
        plot(IPI_sort(keep_mask_sort == 1),inhib_dur_sort(keep_mask_sort == 1))
        hold on
    end
 
%% dot plot showing inhib dur to a single pulse and to the second pulse for different IPIs

    markers = {'.','s'};
    marker_size = [22,8];
    f=figure();
    f.Name = 'Duncan_Han_dblPulse_inhib_vsSinglePulse';
    hold on
    for m = 1:2 % duncan, then han
        keep_mask = monkey_idx == m-1; % matlab indexing causes issues...
        
        plot(inhib_dur(keep_mask,single_pulse_idx),inhib_dur(keep_mask,input_data.num_pulses == 2 & input_data.IPI == 200),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,0),'markerfacecolor',getColorFromList(1,0),'markersize',marker_size(m));
        
        plot(inhib_dur(keep_mask,single_pulse_idx),inhib_dur(keep_mask,input_data.num_pulses == 2 & input_data.IPI == 20),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,1),'markerfacecolor',getColorFromList(1,1),'markersize',marker_size(m));
        
        plot(inhib_dur(keep_mask,single_pulse_idx),inhib_dur(keep_mask,input_data.num_pulses == 2 & input_data.IPI == 10),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,2),'markerfacecolor',getColorFromList(1,2),'markersize',marker_size(m));
        
    end
    
    unity_line = plot([-20,200],[-20,200],'k--','linewidth',1.5);
    double_line = plot([-20,200],[-40,400],'--','color',[0.5,0.5,0.5],'linewidth',1.5);
    
    xlabel('Inhibition duration to single pulse (ms)');
    ylabel('Inhibition duration to second pulse (ms)');
    l=legend('200ms','20ms','10ms'); % actual IPI is [10.6,20.6,201.03]
    set(l,'box','off','fontsize',14,'location','southeast');
    uistack(unity_line,'bottom');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    xlim([0,200]); ylim([0,200])
%     ylim([0,1]);


%% dot plot showing inhib dur due to superposition vs. second pulse inhib dur
  
    markers = {'.','s'};
    marker_size = [22,8];
    f=figure();
    f.Name = 'Duncan_Han_dblPulse_inhib_vsSuperposition';
    hold on
    for m = 1:2 % duncan, then han
        keep_mask = monkey_idx == m-1; % matlab indexing causes issues...
        inhib_dur_superposition = repmat(inhib_dur(:,single_pulse_idx),1,size(inhib_dur,2)) + max(repmat(inhib_dur(:,single_pulse_idx),1,size(inhib_dur,2)) - repmat(input_data.IPI,size(inhib_dur,1),1),zeros(size(inhib_dur)));
        
        plot(inhib_dur_superposition(keep_mask,input_data.num_pulses == 2 & input_data.IPI == 200),inhib_dur(keep_mask,input_data.num_pulses == 2 & input_data.IPI == 200),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,0),'markerfacecolor',getColorFromList(1,0),'markersize',marker_size(m));
        
        plot(inhib_dur_superposition(keep_mask,input_data.num_pulses == 2 & input_data.IPI == 20),inhib_dur(keep_mask,input_data.num_pulses == 2 & input_data.IPI == 20),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,1),'markerfacecolor',getColorFromList(1,1),'markersize',marker_size(m));
        
        plot(inhib_dur_superposition(keep_mask,input_data.num_pulses == 2 & input_data.IPI == 10),inhib_dur(keep_mask,input_data.num_pulses == 2 & input_data.IPI == 10),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,2),'markerfacecolor',getColorFromList(1,2),'markersize',marker_size(m));
        
    end
    
    unity_line = plot([-20,300],[-20,300],'k--','linewidth',1.5);

    xlabel('Inhibition duration to superposition (ms)');
    ylabel('Inhibition duration to second pulse (ms)');
    l=legend('200ms','20ms','10ms'); % actual IPI is [10.6,20.6,201.03]
    set(l,'box','off','fontsize',14,'location','northwest');
    uistack(unity_line,'bottom');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    xlim([0,300]); ylim([0,300])
%     ylim([0,1]);


%% histogram paired difference between single pulse inhib duration and 
    
    single_pulse_idx = find(input_data.num_pulses == 1);
    dbl_pulse_idx = find(input_data.num_pulses == 2);
    bin_edges = [-50:10:200];
    
    figure();
    inhib_dur = [];
    % get relevant inhib dur
    for u = 1:numel(array_data)
%         if(sum(isnan(inhibStruct{u}.inhib_dur([single_pulse_idx,dbl_pulse_idx])'))==0)
            inhib_dur(end+1,:) = inhibStruct{u}.inhib_dur;
%         end
    end
    
    for d = 1:numel(dbl_pulse_idx)
        subplot(numel(dbl_pulse_idx),1,d);
    
        inhib_dur_diff = diff(inhib_dur(:,[1,dbl_pulse_idx(d)])')';
        histogram(inhib_dur_diff,bin_edges);
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

