%% set file names 

    input_data.folderpath = 'E:\Data\Joseph\Han_stim_data\Han_20191112_longTrains_dukeGen2\';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    % input_data.mapFileName = 'mapFileR:\limblab-archive\Retired Animal Logs\Monkeys\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
%     input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';


    input_data.IPI = 1000./[125,100,75,50]; % Hz, -1 means single pulse
    input_data.num_pulses = [125,100,75,50]*4;

    
    input_data.num_conditionitions = numel(input_data.IPI);

    input_data.nom_freq = 2;
    input_data.window = [-500,15000]; % ms before and ms after 1st pulse
    input_data.bin_size = 2; % in ms
    input_data.chan_rec = 8;

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
    
    for fileNumber = 1%:numel(fileList)
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
                array_data{u}.spikeTrialTimes = cell(input_data.num_conditionitions,1);
                array_data{u}.trial_num = cell(input_data.num_conditionitions,1);
                array_data{u}.num_stims = zeros(input_data.num_conditionitions,1);
                array_data{u}.binCounts = cell(input_data.num_conditionitions,1);
                array_data{u}.binEdges = cell(input_data.num_conditionitions,1);
                array_data{u}.kin = cell(input_data.num_conditionitions,1);
                array_data{u}.force = cell(input_data.num_conditionitions,1);
                for condition = 1:numel(input_data.IPI)
                    array_data{u}.kin{condition}.x = [];
                    array_data{u}.kin{condition}.y = [];
                    array_data{u}.kin{condition}.vx = [];
                    array_data{u}.kin{condition}.vy = [];
                    array_data{u}.kin{condition}.ax = [];
                    array_data{u}.kin{condition}.ay = [];
                    array_data{u}.kin{condition}.speed = [];
                    array_data{u}.kin{condition}.mean_speed = [];
                end
            end
        end
        
        % make stimOn a timestamp
        t = (0:1:stimInfo.stimOn(end))/30000;
        stimInfo.stimOn = t(stimInfo.stimOn);

        is_1st_pulse_mask = ones(numel(stimInfo.stimOn),1);
        for st = 1:numel(stimInfo.stimOn)
            if(is_1st_pulse_mask(st) == 1) % look for conditionition based on IPI to next stimulation
                if(st == numel(stimInfo.stimOn)) % last stim was a single pulse
                    conditionition = find(input_data.IPI < 0);
                else
                    IPI = (stimInfo.stimOn(st+1) - stimInfo.stimOn(st))*1000; % in ms
                    num_pulses = find(diff(stimInfo.stimOn(st:end)) > max(input_data.IPI/1000)+0.02,1,'first');
                    if(isempty(num_pulses)) % last pulse
                        num_pulses = numel(stimInfo.stimOn)-st+1;
                    end
                    conditionition_list = find(input_data.num_pulses == num_pulses);
                    [min_diff,conditionition_idx] = min(abs(IPI - input_data.IPI(conditionition_list)));
                    conditionition = conditionition_list(conditionition_idx);
                    if(min_diff > 50)
                        conditionition = find(input_data.IPI < 0);
                    end
                    if(input_data.IPI(conditionition) > 0) % double pulse case
                        is_1st_pulse_mask(st+1:st+input_data.num_pulses(conditionition)-1) = 0;
                    end
                end

                % we know which conditionition and have dealt with the double pulse
                % stuff, so now all we need to do is extract relevant spikes based
                % on input_data.window and store them into array_data

                for u = 1:numel(unit_idx)
                    direct_spike_mask = cds.units(unit_idx(u)).spikes.ts > stimInfo.stimOn(st) + input_data.window(1)/1000 & cds.units(unit_idx(u)).spikes.ts < stimInfo.stimOn(st) + input_data.window(2)/1000;
                    spike_ts = cds.units(unit_idx(u)).spikes.ts(direct_spike_mask == 1) - stimInfo.stimOn(st);
                    array_data{u}.spikeTrialTimes{conditionition} = [array_data{u}.spikeTrialTimes{conditionition}, spike_ts'];
                    array_data{u}.num_stims(conditionition) = array_data{u}.num_stims(conditionition) + 1;
                    array_data{u}.trial_num{conditionition} = [array_data{u}.trial_num{conditionition}, ones(1,numel(spike_ts))*array_data{u}.num_stims(conditionition)];
                    
                    if(isprop(cds,'kin') && ~isempty(cds.kin))
                        kin_win = [find(cds.kin.t > stimInfo.stimOn(st)+input_data.window(1)/1000,1,'first'),...
                            find(cds.kin.t > stimInfo.stimOn(st)+input_data.window(2)/1000,1,'first')];
                        while(diff(kin_win) ~= diff(input_data.window))
                            if(diff(kin_win) < diff(input_data.window))
                                kin_win(end) = kin_win(end)+1;
                            else
                                kin_win(end) = kin_win(end)-1;
                            end
                        end
                        if(numel(kin_win) == 2)
                            array_data{u}.kin{conditionition}.x(end+1,:) = cds.kin.x(kin_win(1):kin_win(2));
                            array_data{u}.kin{conditionition}.y(end+1,:) = cds.kin.y(kin_win(1):kin_win(2));
                            array_data{u}.kin{conditionition}.vx(end+1,:) = cds.kin.vx(kin_win(1):kin_win(2));
                            array_data{u}.kin{conditionition}.vy(end+1,:) = cds.kin.vy(kin_win(1):kin_win(2));
                            array_data{u}.kin{conditionition}.ax(end+1,:) = cds.kin.ax(kin_win(1):kin_win(2));
                            array_data{u}.kin{conditionition}.ay(end+1,:) = cds.kin.ay(kin_win(1):kin_win(2));

                            array_data{u}.kin{conditionition}.speed(end+1,:) = sqrt(cds.kin.vx(kin_win(1):kin_win(2)).^2 + cds.kin.vy(kin_win(1):kin_win(2)).^2);
                        else
                            array_data{u}.kin{conditionition}.x(end+1,:) = nan;
                            array_data{u}.kin{conditionition}.y(end+1,:) = nan;
                            array_data{u}.kin{conditionition}.vx(end+1,:) = nan;
                            array_data{u}.kin{conditionition}.vy(end+1,:) = nan;
                            array_data{u}.kin{conditionition}.ax(end+1,:) = nan;
                            array_data{u}.kin{conditionition}.ay(end+1,:) = nan;

                            array_data{u}.kin{conditionition}.speed(end+1,:) = nan;
                        end
                        
                    end
                    
                    if(isprop(cds,'force') && ~isempty(cds.force))
                        
                        
                    end
                    
                end
            end

        end
        
        % get bin counts for each conditionition
        bin_edges = input_data.window(1):input_data.bin_size:input_data.window(2);

        for u = 1:numel(array_data)
            for condition = 1:numel(array_data{u}.binCounts)
                array_data{u}.binEdges{condition} = bin_edges;
                array_data{u}.binCounts{condition} = histcounts(array_data{u}.spikeTrialTimes{condition}*1000,bin_edges);
                array_data{u}.kin{condition}.mean_speed = mean(array_data{u}.kin{condition}.speed,2);
            end
            
            array_data{u}.monkey = input_data.monkey(7:end);
            array_data{u}.CHAN_LIST = (input_data.chan_rec);
            array_data{u}.stimData = array_data{u}.trial_num;
            array_data{u}.numStims = array_data{u}.num_stims;
        end
        

        
    end

    
    
    
%% rebin 
    binSize = 1; % ms
    
    for u = 1:numel(array_data)
        input_data_all{u}.bin_size = binSize;
        array_data = rebinArrayData(array_data,binSize);
    end

    
%% plot raster 

    for arrIdx = 8%:numel(array_data)   
        array_data{arrIdx}.monkey = 'Han';
        optsPlotFunc.BIN_SIZE = mode(diff(array_data{arrIdx}.binEdges{1,1}));
        optsPlotFunc.FIGURE_SAVE = 0;
        optsPlotFunc.FIGURE_DIR = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\';
        optsPlotFunc.FIGURE_PREFIX = [array_data{arrIdx}.monkey,'_DblPulseTrains_'];

        optsPlotFunc.PRE_TIME = 10/1000;
        optsPlotFunc.POST_TIME = 20/1000;
        optsPlotFunc.SORT_DATA = '';

        optsPlotFunc.STIMULATION_LENGTH = 0.453;

        optsPlotFunc.MARKER_STYLE  = '.';

        optsPlotFunc.PLOT_AFTER_STIMULATION_END = 1;

        rasterPlots = plotRasterStim(array_data{arrIdx},arrIdx,optsPlotFunc);
    end
        

%% plot PSTH for each condition
    for u = 1:numel(array_data)
        input_data.amp_freq = 1;
%         array_data{u}.monkey = 'Han';
%         array_data{u}.num_stims = array_data{u}.numStims;
%       input_data_all{u}.window = [min(array_data{unit_idx}.binEdges{1}),max(array_data{unit_idx}.binEdges{1})];
        input_data.window = [-1500,15000];
        input_data.unit_idx = u;
%         input_data.chan_rec = num2str(array_data{u}.CHAN_LIST);%{1};
        input_data.num_cols = 3;
        input_data.account_for_artifact = 0;
        plotPSTHArrayData(array_data,input_data);
%         f=figure(1);
%         saveFiguresLIB(f,fpath,f.Name);
%         close all
    end
    
    
%% plot PSTH for each condition for low and high speed cases
    mean_speed_cutoff = [0,4;12,10000000]; % plot speeds in range
    input_data.suffix = 'speed0-2';
    input_data.window = [-75,250];
    
    for speed_idx = 1:size(mean_speed_cutoff,1)
        array_data_trim = trimSpeedArrayData(array_data,mean_speed_cutoff(speed_idx,:));
        
        input_data.unit_idx = 1%:1:numel(array_data);
        
        plotPSTHArrayData(array_data_trim,input_data);
    end
    
%% look at response to the train as a function of the baseline activity or speed
    stim_firing_rate = nan(numel(array_data),1);
    baseline_firing_rate = nan(numel(array_data),1);
    
    spike_window = [0,200]; % in ms
    baseline_window = [-80,-5]; % in ms
    
    freq_plot = 20;
    freq_idx = [];
    for unit = [1:numel(array_data)]
%         freq_idx = find(abs(1000/freq_plot - input_data_all{unit}.IPI) < 2 & input_data_all{unit}.num_pulses > 2);
%         freq_idx = freq_idx(1);
        
        baseline_firing_rate(unit) = mean(getSpikesInWindow(array_data{unit},freq_idx,baseline_window,1),'omitnan');
%         stim_firing_rate(unit) = getSpikesInWindow(array_data{unit},freq_idx,spike_window,1);
        
    end

%     plot(baseline_firing_rate,stim_firing_rate,'.','markersize',20)
    histogram(baseline_firing_rate,20)
%% plot response to train against stimulation frequency, use same number of pulses throughout
    markers = {'.','s'};
    marker_size = [22,8];
    
    post_window = [0,5]; % in ms
    num_pulses = 11;
    use_time_for_fit = 1;
    baseline_window = [-80,-5]; % in ms
    slope_window = [10,190]; % in ms
    slope_bin_size = 5;
    
    IPI_plot = [50,20,10,5];

    slopes = nan(numel(array_data),4);
    evoked_response = nan(numel(array_data),4);
    monkey_idx = zeros(numel(array_data),1); % 1 = han, 0 = duncan
    
    num_plot = 10; plot_counter = 0;
    
    for u = 1:numel(array_data)  

        num_pulses_use = input_data_all{u}.num_pulses; % number of pulses to use, deals with situation when we want more pulses than were sent for a conditionition
        num_pulses_use(num_pulses_use > num_pulses) = num_pulses;
                        
        % get baseline firing rate
        baseline_num_spikes = mean(getSpikesInWindow(array_data{u},[],baseline_window,0),'omitnan'); % baseline to subtract is 0, output in firing rate (1)
        
        % define arrays to store evoked data
        stim_firing_rate = zeros(1,numel(array_data{u}.binCounts));
        evoked_per_pulse = zeros(1,numel(array_data{u}.binCounts));
        slopes_temp = nan(1,numel(array_data{u}.binCounts));
        
        for condition = 1:numel(array_data{u}.binCounts)
            if(~isempty(array_data{u}.PULSE_TIMES{condition}))
                per_pulse_response = zeros(num_pulses_use(condition),1);
                time_data = array_data{u}.PULSE_TIMES{condition}{1}*1000;
                for pulse = 1:num_pulses_use(condition)
                    % convert window into indices
                    window_offset = time_data(pulse);
                    
                    per_pulse_response(pulse) = getSpikesInWindow(array_data{u},condition,post_window+window_offset,1) - baseline_num_spikes/diff(baseline_window)*1000; % num_spikes above baseline
                    
                end
                
                if(use_time_for_fit)
                    b_pulse = [ones(numel(per_pulse_response),1),time_data(1:numel(per_pulse_response))']\per_pulse_response;
                else
                    b_pulse = [ones(numel(per_pulse_response),1),(1:1:numel(per_pulse_response))']\per_pulse_response;
                end
                slopes_temp(condition) = b_pulse(2);
                
                % get FR across whole train 
                stim_firing_rate(condition) = sum(per_pulse_response)/(diff(post_window)*num_pulses_use(condition))/1000; 
                
                % get evoked per stim
                evoked_per_pulse(condition) = (stim_firing_rate(condition))/num_pulses_use(condition);
                
            end
        end
        
        % get slope during train
%         b_train = getSlopeDuringTrain(array_data{u},[],slope_window,slope_bin_size);
%         slopes_temp = b_train(2,:);
        
        keep_mask = (input_data_all{u}.num_pulses > 3) & ~isnan(evoked_per_pulse); 
        
        % remove very long train conditionition
        keep_mask(find(input_data_all{u}.num_pulses == 41 & abs(input_data_all{u}.IPI-50) < 1)) = 0;
        
        for i = 1:numel(IPI_plot)
            IPI_idx = find(abs(input_data_all{u}.IPI - IPI_plot(i)) < 1 & input_data_all{u}.num_pulses > 3 & keep_mask == 1);
            if(~isempty(IPI_idx))
                evoked_response(u,i) = evoked_per_pulse(IPI_idx);
                slopes(u,i) = slopes_temp(IPI_idx);
            end
        end
        monkey_idx(u) = strcmpi(array_data{u}.monkey,'Han'); % 1 = han, 0 = duncan
        
%         if(sum(keep_mask) >= 3 && plot_counter < num_plot)
%             plot_counter = plot_counter + 1;
%             % deal with freq_stim not being sorted
%             freq_stim = 1000./input_data_all{u}.IPI; 
%             freq_stim(freq_stim < 0) = 0;
%             [freq_stim_plot,sort_idx] = sort(freq_stim(keep_mask == 1));
%             evoked_per_stim_plot = evoked_per_pulse(keep_mask == 1); 
%             evoked_per_stim_plot = evoked_per_stim_plot(sort_idx);
%             slopes_plot = slopes_temp(keep_mask == 1); 
%             slopes_plot = slopes_plot(sort_idx);
% 
%             % plot data
%             plot(freq_stim_plot,slopes_plot,'marker',markers{monkey_idx(u) + 1},...
%                 'linewidth',1.5,'markersize',marker_size(monkey_idx(u) + 1),'color',getColorFromList(2,plot_counter),'markerfacecolor',getColorFromList(2,plot_counter))
%             
%             xlabel('Stimulation frequency (Hz)')
%             ylabel('\DeltaFR during train (Hz/ms)')
%             formatForLee(gcf)
%             set(gca,'fontsize',14)
%         end
        
    end
    
% plot 20Hz vs higher freq

    f=figure(); hold on;
    f.Name = 'Duncan_Han_following_summary_trains';
    plot([-100,100],[-100,100],'k--','linewidth',1.5)
    hold on
    for idx = 1:2 % duncan, then han
        keep_mask = monkey_idx == idx-1; % matlab indexing causes issues...
        
        plot(slopes(keep_mask,1),slopes(keep_mask,2),'linestyle','none',...
            'marker',markers{idx},'color',getColorFromList(5,2),'markerfacecolor',getColorFromList(5,2),'markersize',marker_size(idx));
%         
        plot(slopes(keep_mask,1),slopes(keep_mask,3),'linestyle','none',...
            'marker',markers{idx},'color',getColorFromList(5,3),'markerfacecolor',getColorFromList(5,3),'markersize',marker_size(idx));
        
        plot(slopes(keep_mask,1),slopes(keep_mask,4),'linestyle','none',...
            'marker',markers{idx},'color',getColorFromList(5,4),'markerfacecolor',getColorFromList(5,4),'markersize',marker_size(idx));
        
    end
    
    if(use_time_for_fit)
        xlabel('Slope during 20Hz trains (Hz/ms)');
        ylabel('Slope during higher freq trains (Hz/ms)');
    else
        xlabel('Slope during 20Hz trains (Hz/pulse)');
        ylabel('Slope during higher freq trains (Hz/pulse)');
    end
%     set(l,'box','off','fontsize',14,'location','northwest');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    xlim([min(slopes(:,1))*1.3,max(slopes(:,1))*1.3]); ylim([-0.5,2])
    ylim([min([slopes(:,2);slopes(:,3);slopes(:,4)])*1.1,max([slopes(:,2);slopes(:,3);slopes(:,4)])*1.1]);
    
%     l=legend('50Hz','94Hz','180Hz')
%     set(l,'box','off')
    

%% stats to see if slopes are more negative than the 20Hz slope
    x_data = []; y_data = []; freqs = [20,50,95,180];
    for j = 1:4
        slopes_use = slopes(~isnan(slopes(:,j)),j);
        x_data = [x_data; freqs(j)*ones(numel(slopes_use),1)];
        y_data = [y_data; slopes_use];
    end
    
    mdl = fitlm(x_data,y_data)
    
    p_list = zeros(3,1);
    for j = 2:4
        p_list(j-1) = ranksum(slopes(:,1),slopes(:,j),'tail','right');
    end
    p_list
    
    
%% plot response to each pulse in the trains of each frequency for each neuron
% 1 plot per frequency
    post_window = [0,5]; % in ms
    max_num_pulses = 41;
    
    baseline_window = [-150,-5]; % in ms
    
    freq_plot = [200,100,50,20];
    
    response_to_each_pulse = nan(numel(array_data),numel(freq_plot),max_num_pulses); % 41 pulses in '200'Hz train
    
    for u = 1:numel(array_data) 
        
        % get baseline firing rate
        baseline_firing_rate = mean(getSpikesInWindow(array_data{u},[],baseline_window,1),'omitnan'); % want in Hz
        
        for condition = 1:numel(array_data{u}.binCounts)
            if(input_data_all{u}.num_pulses(condition) > 2)
                freq_idx = find(abs(input_data_all{u}.IPI(condition)-1000./freq_plot) < 2); 
                    
                for pulse = 1:input_data_all{u}.num_pulses(condition)
                    % convert window into indices
                    window_offset = input_data_all{u}.IPI(condition)*(pulse-1);

                    % store fr to that pulse
                    response_to_each_pulse(u,freq_idx,pulse) = getSpikesInWindow(array_data{u},condition,post_window+window_offset,1)- baseline_firing_rate*diff(post_window)/1000; % count number of spikes above baseline
                    
                end
                
            end
        end
        figure();
        plot(squeeze(response_to_each_pulse(u,:,:))');
    end
        
  
%% compute baseline firing rate for all neurons across all conditions....

    
%% look at inhibition
    optsInhibPlot = [];
    optsInhibPlot.PRE_WINDOW = [-80,-5]; % ms
    post_window = [0,250]; % ms
    optsInhibPlot.MAX_TIME_START = 225; % ms
    optsInhibPlot.BIN_SIZE = 1;
    optsInhibPlot.KERNEL_LENGTH = 10;
    blank_time = [0,5]; % ms
    
    optsInhibPlot.PW1 = 200;
    optsInhibPlot.PW2 = 200;
    optsInhibPlot.POL = 0; % 0 is cathodic first
    optsInhibPlot.NUM_CONSECUTIVE_BINS = 10;
    
    inhibStruct = {};
    
    inhib_dur = [];
    monkey_idx = [];
    for u = 1:numel(array_data)
        if(numel(array_data{u}.num_stims) == 8)
            for cond = 1:numel(array_data{u}.spikeTrialTimes)
                if(~isempty(array_data{u}.PULSE_TIMES{cond}) && numel(array_data{u}.PULSE_TIMES{cond}{1}) == 1) % single pulse
                    optsInhibPlot.POST_WINDOW(cond,:) = post_window;
                    optsInhibPlot.BLANK_TIME(cond,:) = blank_time;
                elseif(~isempty(array_data{u}.PULSE_TIMES{cond}) && numel(array_data{u}.PULSE_TIMES{cond}{1}) == 2)% double pulse
                    optsInhibPlot.POST_WINDOW(cond,:) = post_window + array_data{u}.PULSE_TIMES{cond}{1}(2)*1000;
                    optsInhibPlot.BLANK_TIME(cond,:) = blank_time + array_data{u}.PULSE_TIMES{cond}{1}(2)*1000;
                else % train
                    optsInhibPlot.POST_WINDOW(cond,:) = post_window + 200;
                    optsInhibPlot.BLANK_TIME(cond,:) = blank_time + 200;

                end

            end

            [inhibStruct{u},figure_handles] = plotInhibitionDuration(array_data{u},optsInhibPlot);
            
            if(~isempty(inhibStruct{u}))
                monkey_idx(end+1,1) = strcmpi(array_data{u}.monkey,'Han'); % 1 = han, 0 = duncan 
                inhib_dur(end+1,:) = inhibStruct{u}.inhib_dur;
            end
        end
        % plot inhib duration for each condition, have to sort IPI first
%         [IPI_sort,sort_idx] = sort(input_data.IPI);
%         inhib_dur_sort = inhibStruct{u}.inhib_dur(sort_idx);
%         keep_mask_sort = keep_mask(sort_idx);
%         plot(IPI_sort(keep_mask_sort == 1),inhib_dur_sort(keep_mask_sort == 1))
%         hold on
    end
 
% dot plot showing inhib dur to a single pulse and to the last pulse for
% different train frequencies
    markers = {'.','s'};
    marker_size = [22,8];
    f=figure();
    f.Name = 'Duncan_Han_dblPulse_inhib_trains';
    hold on
    for m = 1:2 % duncan, then han
        keep_mask = monkey_idx == m-1; % matlab indexing causes issues...
        
        plot(inhib_dur(keep_mask,1),inhib_dur(keep_mask,5),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(5,1),'markerfacecolor',getColorFromList(5,1),'markersize',marker_size(m));
        
        plot(inhib_dur(keep_mask,1),inhib_dur(keep_mask,4),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(5,2),'markerfacecolor',getColorFromList(5,2),'markersize',marker_size(m));
        
        plot(inhib_dur(keep_mask,1),inhib_dur(keep_mask,3),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(5,3),'markerfacecolor',getColorFromList(5,3),'markersize',marker_size(m));
        
        plot(inhib_dur(keep_mask,1),inhib_dur(keep_mask,2),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(5,4),'markerfacecolor',getColorFromList(5,4),'markersize',marker_size(m));
        
    end
    
    unity_line = plot([-20,100],[-20,100],'k--','linewidth',1.5);
    
    xlabel('Inhibition duration to single pulse (ms)');
    ylabel('Inhibition duration to last pulse (ms)');
%     l=legend('200ms','20ms','10ms'); % actual IPI is [10.6,20.6,201.03]
%     set(l,'box','off','fontsize',14,'location','southeast');
    uistack(unity_line,'bottom');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    xlim([0,100]); ylim([0,100])
%     ylim([0,1]);

%% plot showing percentage of cells displaying inhibition
    freqs_plot = [-1,20,50,100,180]; % -1 means single pulse
    num_units_total = zeros(numel(freqs_plot),1);
    num_units_inhib = zeros(numel(freqs_plot),1);
    
    for iUnit = 1:numel(inhibStruct)
        if(~isempty(inhibStruct{iUnit}))
            for condition = 1:numel(array_data{iUnit}.spikeTrialTimes)
                if(~isempty(array_data{iUnit}.PULSE_TIMES{condition}))
                    num_stims = numel(array_data{iUnit}.PULSE_TIMES{condition}{1});
                    if(num_stims == 1)
                        idx = 1;
                    elseif(num_stims == 2)
                        idx = [];
                        freq_use = -1; % won't go in next if statement
                    else % train
                        freq_use = 1/(array_data{iUnit}.PULSE_TIMES{condition}{1}(2) - array_data{iUnit}.PULSE_TIMES{condition}{1}(1));
                    end
                    
                    if(num_stims == 1)
                        idx = 1;
                    elseif(num_stims == 2) % do nothing
                        idx = [];
                    elseif(freq_use > 150)
                        idx = 5;
                    elseif(freq_use > 80)
                        idx = 4;
                    elseif(freq_use > 35)
                        idx = 3;
                    elseif(freq_use > 10)
                        idx = 2;
                    end
                    
                    if(~isempty(idx))
                        num_units_total(idx) = num_units_total(idx) + 1;
                        if(inhibStruct{iUnit}.is_inhib(condition))
                            num_units_inhib(idx) = num_units_inhib(idx) + 1;
                        end
                    end
                    
                    
                end
            end
        end
        
    end

    figure();
    hold on
    
    plot(freqs_plot, num_units_inhib./num_units_total);

    
    
%% get decay time constant for all conditions and plot
    is_intermittent = 1;
    bin_size = 50; % ms
    min_rate = 0.5; % Hz
    
    response_amp_time = 500; % ms
    response_amp_num_pulses = -1; % set as a positive number to override time
    response_amp_pulse_window = [1.5,7]; % if using num_pulses, this determines when after each pulse to count spikes
       
    decay_rates = zeros(numel(array_data),12);
    is_responsive = zeros(numel(array_data),12);
    baseline_counts = zeros(numel(array_data),12);
    distance_from_stim = zeros(numel(array_data),1);
    response_amp = zeros(numel(array_data),12);
    response_to_each_pulse = cell(numel(array_data),1);
    chan_stim = []; chan_rec = []; monkey = []; % 1 = han, 0 = duncan
    
    for u = 1:numel(array_data)
        if(strcmpi(array_data{u}.monkey, 'Han')==1)
            map_data = loadMapFile('R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp');
        elseif(strcmpi(array_data{u}.monkey, 'Duncan')==1)
            map_data = loadMapFile('R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp');
        end
        
        % stim window,baseline window,bin size,min mean baseline firing,
        % map_data, only use times when stimulation is occuring
        % (relevant for intermittent stimulation)
        decay_data = getDecayRate(array_data{u},[0,3900],[-2000,-200],bin_size,min_rate*bin_size/1000,map_data,is_intermittent,response_amp_time,response_amp_num_pulses,response_amp_pulse_window); 
        
        decay_rates(u,:) = decay_data.param_list(:,2)';
        is_responsive(u,:) = decay_data.is_responsive';
        baseline_counts(u,:) = decay_data.baseline_counts';
        distance_from_stim(u) = decay_data.distance_from_stim_chan;
        response_amp(u,:) = decay_data.response_amp';
        chan_stim(u,1) = array_data{u}.CHAN_LIST{1};
        chan_rec(u,1) = array_data{u}.CHAN_REC;
        monkey(u,1) = strcmpi(array_data{u}.monkey,'Han');% 1 = han, 0 = duncan
        response_to_each_pulse{u} = decay_data.response_to_each_pulse;
    end
%% get time constant using response to each pulse instead of time data

    for unit = 1:numel(array_ data)
        for condition = 1:12
            
            x_data_stim = (1:210)*1/51;
            y_data_stim = response_to_each_pulse{unit}{condition}(1:numel(x_data_stim));
            
            [fits,gof] = fit(x_data_stim',y_data_stim,'a*exp(-b*x)','StartPoint',[10,1],'Lower',[0,0]);
            decay_rates(unit,condition) = fits.b;
            
        end
    end
    
%% plot time constant
    amp_freq_data = 0;
    marker_size = 8;
    offset = [0.95,1,1.05];
        
    f=figure(); hold on
    f.Name = 'stim_channel_ampfreq_decay_rate';
    
    if(amp_freq_data)
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
            boxplot_wrapper(x_data(condition),(decay_rates(:,condition)),boxplot_params);
        end
        xlabel('Frequency (Hz)');
        ylabel('Decay rate (1/s)')

        xlim([25,155])
        ax = gca;
        ax.XTick = [51,80,104,131];
        formatForLee(gcf);
        set(gca,'fontsize',14);
        ax.XMinorTick = 'off';
    else
        boxplot_params.linewidth = 1.75;
        boxplot_params.box_width = 1.05;
        boxplot_params.whisker_width = boxplot_params.box_width-0.02;
        boxplot_params.outlier_marker = '.';
        boxplot_params.outlier_marker_size = 12;
        boxplot_params.use_log_x_scale = 1;
    
        x_data = [50,50,50,100,100,100,200,200,200,4000,4000,4000].*[1.1,1,0.9,1.1,1,0.9,1.1,1,0.9,1.1,1,0.9];
        color_idx = [2,1,0,2,1,0,2,1,0,2,1,0];
        f.Position = [550 473 890 420];
        ax1=subplot(1,2,1);
        for condition = 1:9
            boxplot_params.outlier_color = getColorFromList(1,color_idx(condition));
            boxplot_params.median_color = getColorFromList(1,color_idx(condition));
            boxplot_params.box_color = getColorFromList(1,color_idx(condition));
            boxplot_params.whisker_color = getColorFromList(1,color_idx(condition));
            boxplot_wrapper(x_data(condition),decay_rates(:,condition),boxplot_params);
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
            boxplot_wrapper(x_data(condition),decay_rates(:,condition),boxplot_params);
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
    end

    
    
    
    
%% plot response_amp vs. distance and fit 
    f=figure();
    f.Name = 'AmpFreq_response_amplitude_distance';
    make_boxplots = 1;
    
    distance_bin_size = 500; %um, box plots only
    dist_spacing = 120;
    boxplot_params.linewidth = 1.75;
    boxplot_params.box_width = 95;
    boxplot_params.whisker_width = boxplot_params.box_width-0.02;
    boxplot_params.outlier_marker = '.';
    boxplot_params.outlier_marker_size = 12;
    boxplot_params.use_log_x_scale = 0;
    
    if(make_boxplots == 0) % dots for each rec:stim pair
        subplot_map = [10,11,12,7,8,9,4,5,6,1,2,3];
        amp_freq_fit = {};
        a_params = zeros(4,3);
        b_params = zeros(4,3);

        x_pred = [0:10:4500];
        for i = 1:12
            ax(i) = subplot(4,3,subplot_map(i));
            plot(distance_from_stim,response_amp(:,i),'.')
            amp_freq_fit{i} = fit(distance_from_stim,response_amp(:,i),'a*exp(-x/b)','StartPoint',[50,4000],'upper',[500,1E5]);
            hold on
            plot(x_pred,feval(amp_freq_fit{i},[0:10:4500]),'r--','linewidth',2)

            a_params(ceil(subplot_map(i)/3),mod(subplot_map(i)-1,3)+1) = amp_freq_fit{i}.a;
            b_params(ceil(subplot_map(i)/3),mod(subplot_map(i)-1,3)+1) = amp_freq_fit{i}.b;
        end

        linkaxes(ax,'xy');
        ylim([-30,150])
    else % make boxplots
      
        subplot_map = [4,4,4,3,3,3,2,2,2,1,1,1];
%         subplot_map = [2,2,2,nan,nan,nan,1,1,1,nan,nan,nan];
        dist_offset = dist_spacing*repmat([-1,0,1],1,4);
        color_idx = [0,1,2,0,1,2,0,1,2,0,1,2];
        distance_bin_edges = (0:distance_bin_size:4000);
        distance_bin_centers = distance_bin_edges(1:end-1) + distance_bin_size/2;
        [~,~,distance_bin_map] = histcounts(distance_from_stim,distance_bin_edges);
        
        for condition = 1:12%[1,2,3,7,8,9]
            ax_list(subplot_map(condition)) = subplot(1,4,subplot_map(condition))
            hold on
            plot([0,max(distance_bin_edges)],[0,0],'k-')
            
            boxplot_params.outlier_color = getColorFromList(1,color_idx(condition));
            boxplot_params.median_color = getColorFromList(1,color_idx(condition));
            boxplot_params.box_color = getColorFromList(1,color_idx(condition));
            boxplot_params.whisker_color = getColorFromList(1,color_idx(condition));
            
            for distance_idx = 1:numel(distance_bin_centers)
                distance_responses = response_amp(distance_bin_map == distance_idx,condition);
                if(numel(distance_responses) > 0)
                    boxplot_wrapper(distance_bin_centers(distance_idx)+dist_offset(condition),...
                        distance_responses,boxplot_params);
                end
            end
            
            formatForLee(gcf)
            set(gca,'fontsize',14)
            if(subplot_map(condition) == 1)
                xlabel('Distance (\mum)');
                ylabel('FR above baseline (Hz)');
            end
        end
        linkaxes(ax_list,'xy');
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
        distance_data = [distance_data;distance_from_stim];
        response_data = [response_data;response_amp(:,condition)];
        amp_data = [amp_data;zeros(numel(distance_from_stim),1)+amps(condition)];
        freq_data = [freq_data;zeros(numel(distance_from_stim),1)+freqs(condition)];
        monkey_data = [monkey_data; monkey];
        chan_stim_data = [chan_stim_data; chan_stim + 100*(monkey==1)];
        chan_rec_data = [chan_rec_data; chan_rec + 100*(monkey==1)];
    end
    
    response_data(response_data < 0) = 0;
    
    modelspec = 'resp ~ dist+amp+freq+chan_rec*chan_stim';
    amp_freq_mdl = {};
    % only use 1 monkey at a time
    for idx = unique(monkey_data)'
        data_table = table(response_data(monkey_data==idx),distance_data(monkey_data==idx),amp_data(monkey_data==idx),...
            freq_data(monkey_data==idx),categorical(chan_rec_data(monkey_data==idx)),categorical(chan_stim_data(monkey_data==idx)),...
            'VariableNames',{'resp','dist','amp','freq','chan_rec','chan_stim'});

        amp_freq_mdl{end+1} = fitlm(data_table,modelspec)
%         amp_freq_mdl{end+1} = fitglm(data_table,modelspec,'Distribution','poisson')
    end

    
%     for condition = 1:12
%         subplot(4,3,subplot_map(condition))
%         hold on
%         y_pred = predict(amp_freq_mdl,[x_pred',amps(condition)+zeros(numel(x_pred),1),freqs(condition)+zeros(numel(x_pred),1)]);
%         
%         plot(x_pred,y_pred,'k--','linewidth',2)
%     end
    
%% get and plot mean spike latency per pulse over the course of the train across all neurons
% must have actual pulse times in arrayData, IPI and num_pulses here is for plotting only
    IPI = 1./[131,131,131,104,104,104,80,80,80,51,51,51]; % in s
    num_pulses = [125,125,125,100,100,100,75,75,75,50,50,50]*4.2;
    pulse_width = 0.453/1000; % s
    latency_data = getSpikeLatencyPerPulse(array_data,pulse_width);
    
% plot data
    marker_size = 8;
    f=figure(); hold on
    f.Name = 'Nonstim_channel_spike_latency';
    f.Position = [251.4000 -49.4000 704 500*input_data.num_cols];
    hold on
    for i = 1:numel(latency_data.spike_times)
        subplot(4,3,i)
        plot(IPI(i)*(0:1:(num_pulses(i)-1)),1000*cellfun(@mean,latency_data.spike_times{i}),...
            'o','color',getColorFromList(1,mod(i-1,3)),'markersize',marker_size)
        formatForLee(gcf);
        set(gca,'fontsize',14);
        if(i==1)
            xlabel('Time post stim (s)');
            ylabel('Mean spike latency (ms)')
        end
        ylim([0,ceil(IPI(i)*1000)]);
    end
    
%% plot direct or indirect PSTH
    latency_divide = 2/1000;

    f=figure();
    f.Name = 'Nonstim_channel_spike_PSTH';
    f.Position = [251.4000 -49.4000 704 500*input_data.num_cols];
    hold on
    for condition = 1:numel(latency_data.spike_times)
        spike_times_direct = [];
        spike_times_indirect = [];
        for pulse = 1:numel(latency_data.spike_times{condition})
            direct_spike_mask = latency_data.spike_times{condition}{pulse} < latency_divide;
            spike_times_adj = latency_data.spike_times{condition}{pulse}(direct_spike_mask) + (pulse-1)*IPI(condition);
            spike_times_direct = [spike_times_direct, spike_times_adj];
            
            indirect_spike_mask = latency_data.spike_times{condition}{pulse} > latency_divide;
            spike_times_adj = latency_data.spike_times{condition}{pulse}(indirect_spike_mask) + (pulse-1)*IPI(condition);
            spike_times_indirect = [spike_times_indirect, spike_times_adj];
        end
        binEdges = 0:0.2:5;
        binCountsDirect= histcounts(spike_times_direct,binEdges);
        binCountsIndirect = histcounts(spike_times_indirect,binEdges);
        
        subplot(4,3,condition)

        plot(binEdges(1:end-1),binCountsDirect,'color',getColorFromList(1,3),'linewidth',2)
        hold on
        plot(binEdges(1:end-1),binCountsIndirect,'color',getColorFromList(1,4),'linewidth',2)
        formatForLee(gcf);
        set(gca,'fontsize',14);
        if(condition == 1)
            l = legend('Direct','Indirect');
            set(l,'box','off')
            xlabel('Time post stim (s)');
            ylabel('Number of spikes');
        end
    end
    
%% resample data w/ smaller number of stimulations

    num_boot = 100;
    num_stims_use = 8; % really should be less than what I actually used
    
    p_list = zeros(numel(array_data),numel(array_data{1}.binCounts),num_boot);
    num_sig = zeros(numel(array_data),numel(array_data{1}.binCounts));
    
    for unit_idx = 1:numel(array_data)
        f=figure();f.Position = [381 -83 1387 856];
        bootstrapped_data = [];
        bootstrapped_firing_rate = zeros(numel(array_data{unit_idx}.binCounts),num_boot,numel(array_data{unit_idx}.binCounts{1}));

        window_post_stim = [0,5];
        pulse_data = {};
        input_data_all{unit_idx} = input_data;
        array_data{unit_idx}.num_stims = array_data{unit_idx}.numStims;
        array_data{unit_idx}.trial_num = array_data{unit_idx}.stimData;
        
        for condition = 1:numel(array_data{unit_idx}.binCounts) % for each condition
            if(input_data_all{unit_idx}.num_pulses(condition) == 1) % get timing of each pulse
                pulse_time = 0;
            else
                pulse_time = [0,input_data_all{unit_idx}.IPI(condition)*(1:input_data_all{unit_idx}.num_pulses(condition)-1)];
            end
            pulse_data{condition}.count = zeros(num_boot,numel(pulse_time));



            for boot = 1:num_boot % for each bootstrap iteration
                % downsample data

                stim_idx_use = datasample(1:array_data{unit_idx}.num_stims(condition),...
                    min(num_stims_use,array_data{unit_idx}.num_stims(condition)),'Replace',false);
                 
                spike_trial_times = array_data{unit_idx}.spikeTrialTimes{condition}(sum(array_data{unit_idx}.trial_num{condition} == stim_idx_use') > 0);
                
                % rebin data and store
                bootstrapped_firing_rate(condition,boot,:) = histcounts(spike_trial_times*1000,array_data{unit_idx}.binEdges{condition})/...
                    (mode(diff(array_data{unit_idx}.binEdges{condition}))/1000)/num_stims_use;
                
                % get p-value for comparing early and late stim
                trial_counts = zeros(num_stims_use,2);
                for trial = 1:num_stims_use
                    spike_trial_times = array_data{unit_idx}.spikeTrialTimes{condition}(array_data{unit_idx}.trial_num{condition} == stim_idx_use(trial)');
                    trial_counts(trial,:) = [sum(spike_trial_times > 0 & spike_trial_times < 0.5), sum(spike_trial_times > 3 & spike_trial_times < 3.5)];
                end
                
                [~,p_list(unit_idx,condition,boot)] = ttest(trial_counts(:,1),trial_counts(:,2));
                
                for p = 1:numel(pulse_time)
                    bin_idx = [find(array_data{unit_idx}.binEdges{condition} > pulse_time(p)+window_post_stim(1),1,'first'),...
                        find(array_data{unit_idx}.binEdges{condition} > pulse_time(p)+window_post_stim(2),1,'first')];

                    pulse_data{condition}.count(boot,p) = mean(bootstrapped_firing_rate(condition,boot,bin_idx(1):bin_idx(2)));
                end
            end

            subplot(3,4,condition)
            x_data = array_data{unit_idx}.binEdges{condition}(1:end-1) + mode(diff(array_data{unit_idx}.binEdges{condition}))/2;
            plot(x_data,squeeze(bootstrapped_firing_rate(condition,1:5:end,:)),'linewidth',1)

            subplot(3,4,condition+4)
            x_data = array_data{unit_idx}.binEdges{condition}(1:end-1) + mode(diff(array_data{unit_idx}.binEdges{condition}))/2;
            plot(x_data,squeeze(bootstrapped_firing_rate(condition,1,:)),'linewidth',2)
            mean_bin_counts = squeeze(mean(bootstrapped_firing_rate(condition,:,:)));

            sorted_bootstrapped_firing_rate = sort(bootstrapped_firing_rate,2);
            idx_use = [floor(num_boot*0.025),ceil(num_boot*0.975)];
            firing_rate_bound = [squeeze(sorted_bootstrapped_firing_rate(condition,idx_use,:))];
    %         
            plot(x_data,mean_bin_counts) % plot the mean
            hold on
            plot(x_data, firing_rate_bound,'--')

            subplot(3,4,condition+8)
            bin_idx = [find(array_data{unit_idx}.binEdges{condition} >= 0,1,'first'),...
                        find(array_data{unit_idx}.binEdges{condition} >= 1000,1,'first'),...
                        find(array_data{unit_idx}.binEdges{condition} >= 2500,1,'first'),...
                        find(array_data{unit_idx}.binEdges{condition} >= 3500,1,'first')];
                    
            histogram(mean(bootstrapped_firing_rate(condition,:,bin_idx(1):bin_idx(2)),3) -...
                mean(bootstrapped_firing_rate(condition,:,bin_idx(3):bin_idx(4)),3));
            
            
            
        end
    end
       
    num_sig = sum(p_list < 0.05,3)

