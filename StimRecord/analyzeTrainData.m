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
                    spike_mask = cds.units(unit_idx(u)).spikes.ts > stimInfo.stimOn(st) + input_data.window(1)/1000 & cds.units(unit_idx(u)).spikes.ts < stimInfo.stimOn(st) + input_data.window(2)/1000;
                    spike_ts = cds.units(unit_idx(u)).spikes.ts(spike_mask == 1) - stimInfo.stimOn(st);
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
    binSize = 200; % ms
    
    for u = 1:numel(array_data)
        input_data_all{u}.bin_size = binSize;
        array_data = rebinArrayData(array_data,binSize);
    end

    
%% plot raster 
    for arrIdx = 1%:numel(array_data)   

        optsPlotFunc.BIN_SIZE = mode(diff(array_data{arrIdx}.binEdges{1,1}));
        optsPlotFunc.FIGURE_SAVE = 0;
        optsPlotFunc.FIGURE_DIR = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\';
        optsPlotFunc.FIGURE_PREFIX = [array_data{arrIdx}.monkey,'_DblPulseTrains_'];

        optsPlotFunc.PRE_TIME = 75/1000;
        optsPlotFunc.POST_TIME = 6000/1000;
        optsPlotFunc.SORT_DATA = '';

        optsPlotFunc.STIMULATION_LENGTH = 0.453;

        optsPlotFunc.MARKER_STYLE  = '.';

        optsPlotFunc.PLOT_AFTER_STIMULATION_END = 1;

        rasterPlots = plotRasterStim(array_data{arrIdx},arrIdx,optsPlotFunc);
    end
        
%% plot PSTH for each condition
    for u = 1:numel(array_data)
        array_data{u}.monkey = 'Duncan';
        array_data{u}.num_stims = array_data{u}.numStims;
%       input_data_all{u}.window = [min(array_data{unit_idx}.binEdges{1}),max(array_data{unit_idx}.binEdges{1})];
        input_data.window = [-4000,12000];
        input_data.unit_idx = u;
        input_data.chan_rec = array_data{u}.CHAN_LIST{1};
        input_data.num_cols = 3;
        
        plotPSTHArrayData(array_data,input_data);
    end
    
    
%% plot PSTH for each conditionition for low and high speed cases
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
    
    spike_window = [0,200]; % in s
    baseline_window = [-80,-5]; % in s
    
    freq_plot = 20;
    
    for unit = 1:numel(array_data)
        freq_idx = find(abs(1000/freq_plot - input_data_all{unit}.IPI) < 2 & input_data_all{unit}.num_pulses > 2);
        freq_idx = freq_idx(1);
        
        baseline_firing_rate(unit) = mean(getSpikesInWindow(array_data{unit},freq_idx,baseline_window,1),'omitnan');
        stim_firing_rate(unit) = getSpikesInWindow(array_data{unit},freq_idx,spike_window,1);
        
    end

    plot(baseline_firing_rate,stim_firing_rate,'.','markersize',20)
   
%% plot response to train against stimulation frequency, use same number of pulses throughout
    markers = {'.','s'};
    marker_size = [22,8];
    
    post_window = [0,5]; % in ms
    num_pulses = 401;

    baseline_window = [-150,-5]; % in ms
    slope_window = [0,190]; % in ms
    slope_bin_size = 5;
    
    IPI_plot = [50,20,10,5];

    figure();
    hold on
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
        slopes_temp = zeros(1,numel(array_data{u}.binCounts));
        
        for condition = 1:numel(array_data{u}.binCounts)
            if(input_data_all{u}.num_pulses(condition) > 2)
                per_pulse_response = zeros(num_pulses_use(condition),1);
                for pulse = 1:num_pulses_use(condition)
                    % convert window into indices
                    window_offset = input_data_all{u}.IPI(condition)*(pulse-1);
                    
                    per_pulse_response(pulse) = getSpikesInWindow(array_data{u},condition,post_window+window_offset,0) - baseline_num_spikes*diff(post_window)/diff(baseline_window); % num_spikes above baseline
                    
                end
                
                % get FR across whole train 
                stim_firing_rate(condition) = sum(per_pulse_response)/(diff(post_window)*num_pulses_use(condition))/1000; 
                
                % get evoked per stim
                evoked_per_pulse(condition) = (stim_firing_rate(condition))/num_pulses_use(condition);
                
            end
        end
        
        % get slope during train
        b_train = getSlopeDuringTrain(array_data{u},[],slope_window,slope_bin_size);
        slopes_temp = b_train(2,:);
        
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
        
        if(sum(keep_mask) >= 3 && plot_counter < num_plot)
            plot_counter = plot_counter + 1;
            disp(u)
            % deal with freq_stim not being sorted
            freq_stim = 1000./input_data_all{u}.IPI; 
            freq_stim(freq_stim < 0) = 0;
            [freq_stim_plot,sort_idx] = sort(freq_stim(keep_mask == 1));
            evoked_per_stim_plot = evoked_per_pulse(keep_mask == 1); 
            evoked_per_stim_plot = evoked_per_stim_plot(sort_idx);
            slopes_plot = slopes_temp(keep_mask == 1); 
            slopes_plot = slopes_plot(sort_idx);

            % plot data
            plot(freq_stim_plot,slopes_plot,'marker',markers{monkey_idx(u) + 1},...
                'linewidth',1.5,'markersize',marker_size(monkey_idx(u) + 1),'color',getColorFromList(2,plot_counter),'markerfacecolor',getColorFromList(2,plot_counter))
            
            xlabel('Stimulation frequency (Hz)')
            ylabel('\DeltaFR during train (Hz/ms)')
            formatForLee(gcf)
            set(gca,'fontsize',14)
        end
        
    end
    
% plot 20Hz vs higher freq

    f=figure();
    f.Name = 'Duncan_Han_following_summary_trains';
    plot([-1,1],[-1,1],'k--','linewidth',1.5)
    hold on
    for m = 1:2 % duncan, then han
        keep_mask = monkey_idx == m-1; % matlab indexing causes issues...
        
        plot(slopes(keep_mask,1),slopes(keep_mask,2),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,0),'markerfacecolor',getColorFromList(1,0),'markersize',marker_size(m));
%         
        plot(slopes(keep_mask,1),slopes(keep_mask,3),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,1),'markerfacecolor',getColorFromList(1,1),'markersize',marker_size(m));
        
        plot(slopes(keep_mask,1),slopes(keep_mask,4),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,2),'markerfacecolor',getColorFromList(1,2),'markersize',marker_size(m));
        
    end
    xlabel('Response to 20Hz trains');
    ylabel('Response to higher freq trains');
%     set(l,'box','off','fontsize',14,'location','northwest');
    formatForLee(gcf);
    set(gca,'fontsize',14);
%     xlim([-0.5,2]); ylim([-0.5,2])
%     ylim([0,1]);

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

    

    
%% look at rebound excitation
% size, peak time, onset, duration
% size = integrate above baseline for the duration
% peak time = peak if it exists
% duration = time from onset to offset

    rebound_input_data.post_stim_window = [20,200];
    rebound_input_data.baseline_window = [-80,-20];
    rebound_input_data.threshold_mult = 2;
    rebound_input_data.num_bins_above_thresh = 5;
    rebound_input_data.bin_size = 5;
    
    for u = 32%1:numel(array_data)
        rebound_input_data.train_length = input_data_all{u}.IPI.*input_data_all{u}.num_pulses;
        rebound_input_data.train_length(rebound_input_data.train_length < 0) = 1; % single pulse condition
        
        reboundExcitation{u} = getReboundExcitationStats(array_data{u},rebound_input_data);
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

