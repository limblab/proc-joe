%% determine filename and input data
    inputData.folderpath = 'D:\Lab\Data\OffsetRT\';
    
%     inputData.folderpath = 'D:\Lab\Data\ReactionTime\Han_20180427_training\';
    inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\right S1 20180919\SN 6251-001804.cmp';
    
    
    inputData.task='taskRToff';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayLeftS1'; 
    inputData.monkey='monkeyJoe';
    inputData.labnum = 6;
    
    pwd=cd;
    cd(inputData.folderpath)
    fileList = dir('*.nev*');
    cd(pwd)
%% load in cds and make trial data
    td_all = [];
    num_trials = 0;
    for fileNumber = 2%:numel(fileList)
        cds = commonDataStructure();
        cds.file2cds([inputData.folderpath fileList(fileNumber).name],inputData.task,inputData.ranBy,...
            inputData.monkey,inputData.labnum,inputData.array1,inputData.mapFileName);
        cd(pwd);
    
        if(~isempty(cds.trials) && size(cds.trials,1) > 1 && sum(cds.trials.result == 'R' | cds.trials.result == 'F') ~= 0)
            % convert cds to trial data
            params.event_list = {'tgtAngle';'tgtOnTime';'isStimTrial';'stimCode';'stimOnset';'stimOffset';'isAudioTrial';'audioCode';'audioOnset';'audioOffset';'isCatchTrial'};
            params.trial_results = {'R','F','A','I'};
            params.extra_time = [1,2];

            params.include_ts = 0;
            params.exclude_units = [0:255];
            td_temp = parseFileByTrial(cds,params);

            % append trial data into a single struct
            for t = 1:numel(td_temp)
                td_temp(t).trial_id = td_temp(t).trial_id + num_trials;
            end
            num_trials = num_trials + size(cds.trials,1);

            td_all = [td_all,td_temp];
        end
    end
    
    % remove trials where go cue occurs after end of trial
    td_all([td_all.idx_goCueTime] > [td_all.idx_endTime]) = [];
    
%% get movement onset for each trial

    [~, td_reward] = getTDidx(td_all, 'result', 'r');
    td_reward = td_reward(~isnan([td_reward.idx_goCueTime]));
    
    % get movement onset
    params.field_idx = 1;
    params.start_idx_offset = 100;
    params.be_aggressive = 1;
    params.which_field = 'acc';

    params.threshold_acc = 35; % absolute threshold on acceleration, using this instead of threshold_mult
    params.min_s = 1;
    params.pre_move_thresh = 50;
    
    params.use_emg = 0;
    params.emg_idx = 13;
    
    td_reward = getMoveOnset(td_reward,params);
    
    % put move-onset back into td_all
    reward_idx = [td_reward.trial_id];
    td_all_rt = td_all;
    for td_reward_idx = 1:numel(td_reward)
        td_all_idx = find([td_all.trial_id] == td_reward(td_reward_idx).trial_id);
        td_all_rt(td_all_idx).idx_movement_on = td_reward(td_reward_idx).idx_movement_on;
    end

%% plot RT's for each audio cue



