%% set file names 

    input_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Han\Han_20190401_dblPulse_trains\chan60stim\';
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
    input_data.window = [-50,450]; % 20ms before and 400ms after 1st pulse
    input_data.bin_size = 2; % in ms
    input_data.chan_rec = 60;

    folderpath = input_data.folderpath; % rest of code uses folderpath currently...may have switched this, not 100% certain

    input_data.task='taskCObump';
    input_data.ranBy='ranByJoseph'; 
    input_data.array1='arrayLeftS1'; 
    input_data.monkey='monkeyDuncan';
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

%% rebin

    input_data.bin_size = 2.5; % ms
    bin_edges = input_data.window(1):input_data.bin_size:input_data.window(2);
    for u = 1:numel(array_data)
        for cond = 1:numel(array_data{u}.binCounts)
            array_data{u}.binEdges{cond} = bin_edges;
            array_data{u}.binCounts{cond} = histcounts(array_data{u}.spikeTrialTimes{cond}*1000,bin_edges);
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
                '_numPulses',num2str(input_data.num_pulses(cond)),'_unitID',num2str(u)];
            plot(array_data{u}.binEdges{cond}(1:end-1)+mean(diff(array_data{u}.binEdges{cond}))/2,...
                array_data{u}.binCounts{cond}/(input_data.bin_size/1000)/array_data{u}.num_stims(cond))
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

    
    




