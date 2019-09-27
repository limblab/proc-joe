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

    input_data.bin_size = 2.5; % ms
    bin_edges = input_data.window(1):input_data.bin_size:input_data.window(2);
    for u = 1:numel(array_data)
        for cond = 1:numel(array_data{u}.binCounts)
            array_data{u}.binEdges{cond} = bin_edges;
            array_data{u}.binCounts{cond} = histcounts(array_data{u}.spikeTrialTimes{cond}*1000,bin_edges);
        end
    end
    
     
%% plot raster 
    for arrIdx = 3%1:numel(array_data)   

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
    for unit_idx = 8%1:numel(array_data)
        input_data.unit_idx = unit_idx;
        input_data.chan_rec = array_data{unit_idx}.CHAN_LIST;
        plotPSTHArrayData(array_data,input_data);
    end
    
    
%% count spikes in window after each stimulation, plot those values
    
    window = [1,5]; % in ms
    prob_spike = nan(numel(array_data),input_data.num_conditions,2);
    
    keep_mask = input_data.num_pulses == 2;
    
    
    for u = 1:numel(array_data)
        for cond = 1:numel(array_data{u}.spikeTrialTimes)
            if(keep_mask(cond) == 1)
                
                window_adj = window + input_data.IPI(cond);
                                
                prob_spike(u,cond,1) = sum(array_data{u}.spikeTrialTimes{cond}*1000 > window(1) & ...
                    array_data{u}.spikeTrialTimes{cond}*1000 < window(2))/array_data{u}.num_stims(cond);
                
                prob_spike(u,cond,2) = sum(array_data{u}.spikeTrialTimes{cond}*1000 > window_adj(1) & ...
                    array_data{u}.spikeTrialTimes{cond}*1000 < window_adj(2))/array_data{u}.num_stims(cond);
            end
        end
    end
    
    [sort_IPI,sort_idx] = sort(input_data.IPI);
    prob_spike_sort = prob_spike(:,sort_idx,:);
    keep_mask_sort = keep_mask(sort_idx);
    
  
    figure();
    plot(sort_IPI(keep_mask_sort == 1),prob_spike_sort(:,keep_mask_sort == 1,2) - prob_spike_sort(:,keep_mask_sort == 1,1))


%% look at duration of inhibition after second pulse compared to the single pulse case
    optsInhibPlot = [];
    optsInhibPlot.PRE_WINDOW = [-200,-5];
    post_window = [0,200]; % ms
    optsInhibPlot.MAX_TIME_START = 40; % ms
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
    
    for u = 1:numel(array_data)
        for cond = 1:numel(array_data{u}.spikeTrialTimes)
            if(input_data.num_pulses(cond) == 1) % single pulse
                optsInhibPlot.POST_WINDOW = post_window;
                optsInhibPlot.BLANK_TIME = blank_time;
            elseif(input_data.num_pulses(cond) == 2)% double pulse
                optsInhibPlot.POST_WINDOW = post_window + input_data.IPI(cond);
                optsInhibPlot.BLANK_TIME = blank_time;
                optsInhibPlot.BLANK_TIME(2) = blank_time(2) + input_data.IPI(cond);
            else % train
                optsInhibPlot.POST_WINDOW = post_window + 200;
                optsInhibPlot.BLANK_TIME = blank_time;
                optsInhibPlot.BLANK_TIME(2) = blank_time(2) + 200;

            end

            [inhibStruct{u},figure_handles] = plotInhibitionDuration(array_data{u},optsInhibPlot);

        
        end
        
        % plot inhib duration for each condition, have to sort IPI first
        [IPI_sort,sort_idx] = sort(input_data.IPI);
        inhib_dur_sort = inhibStruct{u}.inhib_dur(sort_idx);
        keep_mask_sort = keep_mask(sort_idx);
        plot(IPI_sort(keep_mask_sort == 1),inhib_dur_sort(keep_mask_sort == 1))
        hold on
    end
    
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
    
        inhib_dur_diff = diff(inhib_dur(:,[1,d+1])')';
        histogram(inhib_dur_diff,bin_edges);
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

