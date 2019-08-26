%% set file names 

    input_data.folderpath = 'C:\Users\jts3256\Desktop\Duncan_stim_data\Duncan_20190327_trains\';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    % input_data.mapFileName = 'mapFileR:\limblab-archive\Retired Animal Logs\Monkeys\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
%     input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';


    input_data.IPI = [-1,5,10,20,50,200,10,20]; % Hz, -1 means single pulse
    input_data.num_pulses = [1,41,21,11,5,2,2,2];
%     input_data.IPI = [-1,10,20,50];
%     input_data.num_pulses = [1,21,11,5];
    
    input_data.train_length = 0.2; % s
    input_data.num_conditions = numel(input_data.IPI);

    input_data.nom_freq = 2;
    input_data.window = [-250,500]; % 20ms before and 400ms after 1st pulse
    input_data.bin_size = 2; % in ms
    input_data.chan_rec = 1;

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
                    num_pulses = sum(stimInfo.stimOn(st:end) < stimInfo.stimOn(st)+1/input_data.nom_freq-0.01);
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
        end
    end

    
    
    
%% rebin 
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
    input_data.unit_idx = unit_idx;
    plotPSTHArrayData(array_data,input_data);

    
    
%% plot PSTH for each condition for low and high speed cases
    mean_speed_cutoff = [0,5]; % plot speeds in range
    input_data.suffix = 'speed3-12';
    
    array_data_trim = array_data;
    for u = 1:numel(array_data)
        for cond = 1:numel(array_data{u}.binCounts)
            trial_speeds = array_data{u}.kin{cond}.mean_speed(array_data{u}.trial_num{cond});
            keep_mask = trial_speeds > mean_speed_cutoff(1) & trial_speeds < mean_speed_cutoff(2);
            
            spike_trial_times = array_data{u}.spikeTrialTimes{cond}(keep_mask==1);
            array_data_trim{u}.binCounts{cond} = histcounts(spike_trial_times*1000,array_data_trim{u}.binEdges{cond});
            array_data_trim{u}.num_stims(cond) = sum(array_data{u}.kin{cond}.mean_speed > mean_speed_cutoff(1) & array_data{u}.kin{cond}.mean_speed < mean_speed_cutoff(2));

        end
    end
    input_data.unit_idx = unit_idx;
    plotPSTHArrayData(array_data_trim,input_data);

    
    
    
    
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
       


