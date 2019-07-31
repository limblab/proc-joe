%% set file names 

    input_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Han\Han_20190730_trains_differentIPIs\';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    % input_data.mapFileName = 'mapFileR:\limblab-archive\Retired Animal Logs\Monkeys\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
%     input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';


    input_data.IPI = [500,2000,20000]; % time between trains
%     input_data.IPI = [-1,10,20,50];
%     input_data.num_pulses = [1,21,11,5];
    
    input_data.train_length = 0.2; % s
    input_data.num_conditions = numel(input_data.IPI);

    input_data.train_freq = 200; %Hz
    input_data.window = [-400,800]; % 20ms before and 400ms after 1st pulse
    
    input_data.bin_size = 2; % in ms
    input_data.chan_rec = 21;

    folderpath = input_data.folderpath; % rest of code uses folderpath currently...may have switched this, not 100% certain

    input_data.task='taskBD';
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

        is_1st_pulse_mask = zeros(numel(stimInfo.stimOn),1); 
        is_1st_pulse_mask(1) = 1;
        new_train = 0;
        for st = 2:numel(stimInfo.stimOn)
            IPI = (stimInfo.stimOn(st) - stimInfo.stimOn(st-1))*1000; % in ms
            if(IPI > 1000/input_data.train_freq + 10) % add 10 ms to be sure
                % we found a transition between trains, 
                is_1st_pulse_mask(st) = 1;
            end
            
            [min_diff,condition] = min(abs(IPI - input_data.IPI));


            % we know which condition and have dealt with the double pulse
            % stuff, so now all we need to do is extract relevant spikes based
            % on input_data.window and store them into array_data
            
            % only take spike data if we are at the start of a train
            if(is_1st_pulse_mask(st) == 1)
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

%% rebin and downsample number of stimulations (if desired)
    input_data.bin_size = 5; % ms
    bin_edges = input_data.window(1):input_data.bin_size:input_data.window(2);
    for u = 1:numel(array_data)
        for cond = 1:numel(array_data{u}.binCounts)
            spike_trial_times = array_data{u}.spikeTrialTimes{cond};
            array_data{u}.binEdges{cond} = bin_edges;
            array_data{u}.binCounts{cond} = histcounts(spike_trial_times*1000,bin_edges);
        end
    end

%% plot PSTH for each condition
    max_y_lim = 0;
    for u = 1%:numel(unit_idx)
        max_y_lim = 0;
        f_list = {};
        for cond = 1:numel(array_data{u}.binCounts)
            f_list{cond} = figure();
            f_list{cond}.Name = ['Han_20190401_chan',num2str(input_data.chan_rec),'rec_IPI',num2str(input_data.IPI(cond)),...
                '_unitID',num2str(u)];
            if(exist('downsample_stims') && downsample_stims == 1) % use num_stims_use 
                plot(array_data{u}.binEdges{cond}(1:end-1)+mean(diff(array_data{u}.binEdges{cond}))/2,...
                    array_data{u}.binCounts{cond}/(input_data.bin_size/1000)/min(num_stims_use,array_data{u}.num_stims(cond)))
            else % use num_stims from array_data
                plot(array_data{u}.binEdges{cond}(1:end-1)+mean(diff(array_data{u}.binEdges{cond}))/2,...
                    array_data{u}.binCounts{cond}/(input_data.bin_size/1000)/array_data{u}.num_stims(cond))
            end
            formatForLee(gcf)
            xlabel('Time after 1st pulse (ms)');
            ylabel('Firing rate (spks/s)')
            set(gca,'fontsize',14);
            max_y_lim = max(max_y_lim,f_list{cond}.Children.YLim(end));
        end
        
        for cond = 1:numel(array_data{u}.binCounts)
            f_list{cond}.Children.YLim(end) = max_y_lim;
        end
    end

    
    
%% bootstrap data w/ smaller number of stimulations

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
       


