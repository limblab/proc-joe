%% script to process reaction time data 
%% determine filename and input data
    inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Han_20180501\';
%     inputData.folderpath = 'D:\Lab\Data\ReactionTime\Han_20180427_training\';
    inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

    inputData.task='taskRT';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayLeftS1'; 
    inputData.monkey='monkeyHan';
    inputData.labnum = 6;
    
    pwd=cd;
    cd(inputData.folderpath)
    fileList = dir('*.nev*');
    cd(pwd)
%% load in cds and extract data

    cds = commonDataStructure();
    cds.file2cds([inputData.folderpath fileList(1).name],inputData.task,inputData.ranBy,...
        inputData.monkey,inputData.labnum,inputData.array1,inputData.mapFileName,'recoverPreSync');
    cd(pwd);
    
%% convert cds to trial data
    params.event_list = {'bumpTime';'bumpMagnitude';'bumpDir';'isStimTrial';'stimCode'};
    params.trial_results = {'R','F'};
%     params.include_stim_times = 1;
%     params.include_true_bump_times = 1;
    td_all = parseFileByTrial(cds,params);
    td_all = getGoCueTime(td_all,cds);
    
    [~, td_reward] = getTDidx(td_all, 'result', 'r');
    td_reward = td_reward(~isnan([td_reward.idx_goCueTime]));
    
%% get movement onset
    params.which_field = 'vel';
    params.field_idx = 1;
    params.start_idx_offset = 15;
    params.threshold_mult = 0.2;
    td_reward = getMoveOnsetAndPeak(td_reward,params);
    
% put movement on back into td_all
    counter = 1;
    reward_idx = [td_reward.trial_id];
    for r_idx = reward_idx
        td_all(r_idx).idx_movement_on = td_reward(counter).idx_movement_on;
        counter = counter + 1;
    end
            
%% plot a set of reaches aligned to go cue with reaction time markers
    opts.MAX_PLOT = 4;
    opts.WHICH_FIELD ='acc';
    opts.WHICH_IDX = 1;
    plotReachesTD(td_reward,opts);
    
%% plot reaction times for each cue, psychometric curve
    opts = [];
    opts.FOLDER_PATH = inputData.folderpath;
    opts.SAVE_FIGURES = 1;
    opts.FIGURE_PREFIX = 'Han_20180501'; % no _ required
    [~,plots] = plotReactionTimeDataTD(td_reward,td_all,opts);


    
    
    
    