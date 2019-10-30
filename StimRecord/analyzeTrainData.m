%% set file names 

    input_data.folderpath = 'C:\Users\jts3256\Desktop\Han_stim_data\Han_20191011_dukeProjBox_trains\';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    % input_data.mapFileName = 'mapFileR:\limblab-archive\Retired Animal Logs\Monkeys\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
%     input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';


    input_data.IPI = [-1,5,10,20,50,50]; % Hz, -1 means single pulse
    input_data.num_pulses = [1,41,21,11,5,41];

    
    input_data.num_conditions = numel(input_data.IPI);

    input_data.nom_freq = 2;
    input_data.window = [-250,2500]; % ms before and ms after 1st pulse
    input_data.bin_size = 2; % in ms
    input_data.chan_rec = 51;

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
                array_data{u}.kin = cell(input_data.num_conditions,1);
                array_data{u}.force = cell(input_data.num_conditions,1);
                for cond = 1:numel(input_data.IPI)
                    array_data{u}.kin{cond}.x = [];
                    array_data{u}.kin{cond}.y = [];
                    array_data{u}.kin{cond}.vx = [];
                    array_data{u}.kin{cond}.vy = [];
                    array_data{u}.kin{cond}.ax = [];
                    array_data{u}.kin{cond}.ay = [];
                    array_data{u}.kin{cond}.speed = [];
                    array_data{u}.kin{cond}.mean_speed = [];
                end
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
                    IPI = (stimInfo.stimOn(st+1) - stimInfo.stimOn(st))*1000; % in ms
                    num_pulses = find(diff(stimInfo.stimOn(st:end)) > max(input_data.IPI/1000)+0.02,1,'first');
                    if(isempty(num_pulses)) % last pulse
                        num_pulses = numel(stimInfo.stimOn)-st+1;
                    end
                    condition_list = find(input_data.num_pulses == num_pulses);
                    [min_diff,condition_idx] = min(abs(IPI - input_data.IPI(condition_list)));
                    condition = condition_list(condition_idx);
                    if(min_diff > 50)
                        condition = find(input_data.IPI < 0);
                    end
                    if(input_data.IPI(condition) > 0) % double pulse case
                        is_1st_pulse_mask(st+1:st+input_data.num_pulses(condition)-1) = 0;
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
                            array_data{u}.kin{condition}.x(end+1,:) = cds.kin.x(kin_win(1):kin_win(2));
                            array_data{u}.kin{condition}.y(end+1,:) = cds.kin.y(kin_win(1):kin_win(2));
                            array_data{u}.kin{condition}.vx(end+1,:) = cds.kin.vx(kin_win(1):kin_win(2));
                            array_data{u}.kin{condition}.vy(end+1,:) = cds.kin.vy(kin_win(1):kin_win(2));
                            array_data{u}.kin{condition}.ax(end+1,:) = cds.kin.ax(kin_win(1):kin_win(2));
                            array_data{u}.kin{condition}.ay(end+1,:) = cds.kin.ay(kin_win(1):kin_win(2));

                            array_data{u}.kin{condition}.speed(end+1,:) = sqrt(cds.kin.vx(kin_win(1):kin_win(2)).^2 + cds.kin.vy(kin_win(1):kin_win(2)).^2);
                        else
                            array_data{u}.kin{condition}.x(end+1,:) = nan;
                            array_data{u}.kin{condition}.y(end+1,:) = nan;
                            array_data{u}.kin{condition}.vx(end+1,:) = nan;
                            array_data{u}.kin{condition}.vy(end+1,:) = nan;
                            array_data{u}.kin{condition}.ax(end+1,:) = nan;
                            array_data{u}.kin{condition}.ay(end+1,:) = nan;

                            array_data{u}.kin{condition}.speed(end+1,:) = nan;
                        end
                        
                    end
                    
                    if(isprop(cds,'force') && ~isempty(cds.force))
                        
                        
                    end
                    
                end
            end

        end
        
        % get bin counts for each condition
        bin_edges = input_data.window(1):input_data.bin_size:input_data.window(2);

        for u = 1:numel(array_data)
            for cond = 1:numel(array_data{u}.binCounts)
                array_data{u}.binEdges{cond} = bin_edges;
                array_data{u}.binCounts{cond} = histcounts(array_data{u}.spikeTrialTimes{cond}*1000,bin_edges);
                array_data{u}.kin{cond}.mean_speed = mean(array_data{u}.kin{cond}.speed,2);
            end
            
            array_data{u}.monkey = input_data.monkey(7:end);
            array_data{u}.CHAN_LIST = (input_data.chan_rec);
            array_data{u}.stimData = array_data{u}.trial_num;
            array_data{u}.numStims = array_data{u}.num_stims;
        end
        

        
    end

    
    
    
%% rebin 
    input_data.bin_size = 2.5; % ms
    bin_edges = input_data.window(1):input_data.bin_size:input_data.window(2);
    for u = 1:numel(array_data)
        for cond = 1:numel(array_data{u}.binCounts)
            spike_trial_times = array_data{u}.spikeTrialTimes{cond};
            array_data{u}.binEdges{cond} = bin_edges;
            array_data{u}.binCounts{cond} = histcounts(spike_trial_times*1000,bin_edges);
        end
    end

    
%% plot raster 
    for arrIdx = 1:numel(array_data)   

        optsPlotFunc.BIN_SIZE = mode(diff(array_data{arrIdx}.binEdges{1,1}));
        optsPlotFunc.FIGURE_SAVE = 0;
        optsPlotFunc.FIGURE_DIR = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\';
        optsPlotFunc.FIGURE_PREFIX = [array_data{arrIdx}.monkey,'_DblPulseTrains_'];

        optsPlotFunc.PRE_TIME = 75/1000;
        optsPlotFunc.POST_TIME = 250/1000;
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
    
    
%% plot PSTH for each condition for low and high speed cases
    mean_speed_cutoff = [12,10000000]; % plot speeds in range
    input_data.suffix = 'speed0-2';
    
    array_data_trim = array_data;
    for u = 1:numel(array_data)
        for cond = 1:numel(array_data{u}.binCounts)
            trial_speeds = array_data{u}.mean_speed{cond}(array_data{u}.trial_num{cond});
            keep_mask = trial_speeds > mean_speed_cutoff(1) & trial_speeds < mean_speed_cutoff(2);
            trials_keep = unique(array_data{u}.trial_num{cond}(keep_mask));
            
            
            spike_trial_times = array_data{u}.spikeTrialTimes{cond}(keep_mask==1);
            array_data_trim{u}.binCounts{cond} = histcounts(spike_trial_times*1000,array_data_trim{u}.binEdges{cond});
            array_data_trim{u}.num_stims(cond) = sum(array_data{u}.mean_speed{cond} > mean_speed_cutoff(1) & array_data{u}.mean_speed{cond} < mean_speed_cutoff(2));
            array_data_trim{u}.speed{cond} = array_data{u}.speed{cond}(trials_keep,:);
            array_data_trim{u}.mean_speed{cond} = array_data{u}.mean_speed{cond}(trials_keep,:);
        end
    end
    input_data.unit_idx = 1;
    plotPSTHArrayData(array_data_trim,input_data);

    
%% look at response to the train as a function of the baseline activity or speed
    num_spikes_per_trial = cell(numel(array_data{1}.num_stims),1);
    num_spikes_baseline = cell(numel(array_data{1}.num_stims),1);
    mean_speed_per_trial = cell(numel(array_data{1}.num_stims),1);
    spike_window = [0,200]/1000; % in s
    baseline_window = [-250,-5]/1000; % in s
    
    for cond = [2,3,4,5]%1:numel(array_data{1}.num_stims)
        for tr = 1:array_data{1}.num_stims(cond)
            % get spikes during the train
            spike_mask = array_data{1}.spikeTrialTimes{cond} > spike_window(1) & array_data{1}.spikeTrialTimes{cond} < spike_window(2) & ...
                array_data{1}.trial_num{cond} == tr;
            num_spikes_per_trial{cond}(end+1,1) = sum(spike_mask);
            mean_speed_per_trial{cond}(end+1,1) = array_data{1}.mean_speed{cond}(tr);
            
            % get baseline spike rate
            spike_mask = array_data{1}.spikeTrialTimes{cond} > baseline_window(1) & array_data{1}.spikeTrialTimes{cond} < baseline_window(2) & ...
                array_data{1}.trial_num{cond} == tr;
            num_spikes_baseline{cond}(end+1,1) = sum(spike_mask);  
            
        end
        figure();
        plot(mean_speed_per_trial{cond},num_spikes_per_trial{cond},'.','markersize',12);
        xlabel('speed (cm/s)');
        ylabel('Num spikes evoked');
        formatForLee(gcf);
        set(gca,'fontsize',14)
        
        figure();
        plot(num_spikes_baseline{cond},num_spikes_per_trial{cond},'.','markersize',12);
        xlabel('Num spikes baseline');
        ylabel('Num spikes evoked');
        formatForLee(gcf);
        set(gca,'fontsize',14)
    end

    
    
%% plot response to train against stimulation frequency
    markers = {'.','s'};
    marker_size = [22,8];
    
    post_window = [0,215]; % in ms
    baseline_window = [-150,-5]; % in ms
    
    input_data.IPI = [-1,5.6,10.8,20,50,200,10,20];
    freq_stim = [1,180,90,50,20,5,100,50];
    figure();
    hold on
    slopes = [];
    norm_resp = zeros(numel(array_data),4);
    monkey_idx = zeros(numel(array_data),1); % 1 = han, 0 = duncan
    
    num_plot = 15; plot_counter = 0;
    
    for u = 3:numel(array_data)        
        num_pulses = [];
        baseline_firing_rate = [];
        stim_firing_rate = [];
        evoked_per_stim = [];
        
        
        % convert window into indices
        post_stim_idx = [find(array_data{u}.binEdges{1} >= post_window(1),1,'first'),...
            find(array_data{u}.binEdges{1} >= post_window(2),1,'first')];
        baseline_idx = [find(array_data{u}.binEdges{1} >= baseline_window(1),1,'first'),...
            find(array_data{u}.binEdges{1} >=baseline_window(2),1,'first')];
        
        for cond = 1:numel(array_data{u}.binCounts)
%             freq_stim(cond) = 1000/input_data.IPI(cond);
            num_pulses(cond) = input_data.num_pulses(cond);
            
            if(num_pulses(cond) == 1) % fix single pulse condition
                freq_stim(cond) = 1;
            end
            
            baseline_firing_rate(cond) = sum(array_data{u}.binCounts{cond}(baseline_idx(1):baseline_idx(2)))/...
                array_data{u}.num_stims(cond)/(diff(baseline_window)-num_pulses(cond))*1000; % convert to Hz
            stim_firing_rate(cond) = sum(array_data{u}.binCounts{cond}(post_stim_idx(1):post_stim_idx(2)))/...
                array_data{u}.num_stims(cond)/(diff(post_window)-num_pulses(cond))*1000; % convert to Hz, subtract 1 ms per artifact
            stim_fit(cond,:) = [ones(diff(post_stim_idx)+1,1), array_data{u}.binEdges{cond}(post_stim_idx(1):post_stim_idx(2))']\...
                array_data{u}.binCounts{cond}(post_stim_idx(1):post_stim_idx(2))';
            
            stim_slope(cond) = stim_fit(cond,2);
            evoked_per_stim(cond) = (stim_firing_rate(cond)-baseline_firing_rate(cond))/num_pulses(cond);
        end
        keep_mask = (num_pulses > 3) & ~isnan(evoked_per_stim);
        
        % get 5 Hz, 100Hz, 200Hz response for future plot
        norm_resp(u,:) = evoked_per_stim([find(freq_stim == 20 & num_pulses > 3,1,'first'),...
            find(freq_stim == 50 & num_pulses > 3,1,'first'),...
            find(freq_stim == 90 & num_pulses > 3,1,'first'),...
            find(freq_stim == 180 & num_pulses > 3,1,'first')]);
        monkey_idx(u) = strcmpi(array_data{u}.monkey,'Han'); % 1 = han, 0 = duncan
        
        if(sum(keep_mask) >= 4 && plot_counter < num_plot)
            plot_counter = plot_counter + 1;
            
            % deal with freq_stim not being sorted
            [freq_stim_plot,sort_idx] = sort(freq_stim(keep_mask == 1));
            evoked_per_stim_plot = evoked_per_stim(keep_mask == 1); 
            evoked_per_stim_plot = evoked_per_stim_plot(sort_idx);
            stim_slope_plot = stim_slope(keep_mask == 1);
            stim_slope_plot = stim_slope_plot(sort_idx);

            % plot data
            plot(freq_stim_plot,evoked_per_stim_plot,'marker',markers{monkey_idx(u) + 1},...
                'linewidth',1.5,'markersize',marker_size(monkey_idx(u) + 1),'color',getColorFromList(2,plot_counter),'markerfacecolor',getColorFromList(2,plot_counter))
            
            xlabel('Stimulation frequency (Hz)')
            ylabel('(F_e - F_b)/N')
            formatForLee(gcf)
            set(gca,'fontsize',14)
        end
        
        % fit following summary and get slopes
        slope_unit = [ones(sum(keep_mask),1),freq_stim(keep_mask == 1)']\evoked_per_stim(keep_mask == 1)';
        slopes(u) = slope_unit(2);
    end
    
    
%% plot 20Hz resp vs. 100 and 200Hz   

    f=figure();
    f.Name = 'Duncan_Han_following_summary_trains';
    plot([-2,8],[-2,8],'k--','linewidth',1.5)
    hold on
    for m = 1:2 % duncan, then han
        keep_mask = monkey_idx == m-1; % matlab indexing causes issues...
        
        plot(norm_resp(keep_mask,1),norm_resp(keep_mask,2),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,0),'markerfacecolor',getColorFromList(1,0),'markersize',marker_size(m));
        
        plot(norm_resp(keep_mask,1),norm_resp(keep_mask,3),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,1),'markerfacecolor',getColorFromList(1,1),'markersize',marker_size(m));
        
        plot(norm_resp(keep_mask,1),norm_resp(keep_mask,4),'linestyle','none',...
            'marker',markers{m},'color',getColorFromList(1,2),'markerfacecolor',getColorFromList(1,2),'markersize',marker_size(m));
        
    end
    xlabel('Response to 20Hz trains');
    ylabel('Response to higher freq trains');
    set(l,'box','off','fontsize',14,'location','northwest');
    formatForLee(gcf);
    set(gca,'fontsize',14);
%     xlim([-0.5,2]); ylim([-0.5,2])
%     ylim([0,1]);

    
%% resample data w/ smaller number of stimulations

    unit_idx = 1;
    num_boot = 1000;
    num_stims_use = 35; % really should be less than what I actually used
    
    bootstrapped_data = [];
    bootstrapped_firing_rate = zeros(numel(array_data{unit_idx}.binCounts),num_boot,numel(array_data{unit_idx}.binCounts{1}));
    
    window_post_stim = [1,10];
    pulse_data = {};
    for cond = 1:numel(array_data{unit_idx}.binCounts) % for each condition
        if(input_data.num_pulses(cond) == 1) % get timing of each pulse
            pulse_time = 0;
        else
            pulse_time = [0,input_data.IPI(cond)*(1:input_data.num_pulses(cond)-1)];
        end
        pulse_data{cond}.count = zeros(num_boot,numel(pulse_time));
        
        
        
        for boot = 1:num_boot % for each bootstrap iteration
            % downsample data
            stim_idx_use = datasample(1:array_data{unit_idx}.num_stims(cond),...
                min(num_stims_use,array_data{unit_idx}.num_stims(cond)),'Replace',false);
            spike_trial_times = array_data{unit_idx}.spikeTrialTimes{cond}(sum(array_data{unit_idx}.trial_num{cond} == stim_idx_use') > 0);
            % rebin data and store
            bootstrapped_firing_rate(cond,boot,:) = histcounts(spike_trial_times*1000,array_data{unit_idx}.binEdges{cond})/...
                (input_data.bin_size/1000)/num_stims_use;
            
            for p = 1:numel(pulse_time)
                bin_idx = [find(array_data{unit_idx}.binEdges{cond} > pulse_time(p)+window_post_stim(1),1,'first'),...
                    find(array_data{unit_idx}.binEdges{cond} > pulse_time(p)+window_post_stim(2),1,'first')];
                
                pulse_data{cond}.count(boot,p) = mean(bootstrapped_firing_rate(cond,boot,bin_idx(1):bin_idx(2)));
            end
        end
        
        figure();
        subplot(2,1,1)
        x_data = array_data{unit_idx}.binEdges{cond}(1:end-1) + mode(diff(array_data{unit_idx}.binEdges{cond}))/2;
        plot(x_data,squeeze(bootstrapped_firing_rate(cond,1,:)),'linewidth',2)
%         mean_bin_counts = squeeze(mean(bootstrapped_firing_rate(cond,:,:)));
% 
%         sorted_bootstrapped_firing_rate = sort(bootstrapped_firing_rate,2);
%         idx_use = [floor(num_boot*0.025),ceil(num_boot*0.975)];
%         firing_rate_bound = [squeeze(sorted_bootstrapped_firing_rate(cond,idx_use,:))];
%         
%         plot(x_data,mean_bin_counts) % plot the mean
%         hold on
%         plot(x_data, firing_rate_bound,'--')
%         subplot(2,1,2)
%         histogram(pulse_data{cond}.count(:,1) - pulse_data{cond}.count(:,end));

    end
       


