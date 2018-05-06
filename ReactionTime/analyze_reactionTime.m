%% script to process reaction time data 
%% determine filename and input data
%     inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Han_20180501\';
    inputData.folderpath = 'D:\Lab\Data\ReactionTime\Han_20180501_training\';
    inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

    inputData.task='taskRT';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayLeftS1'; 
    inputData.monkey='monkeyHan';
    inputData.labnum = 6;
    
    pwd=cd;
    cd(inputData.folderpath)
    fileList = dir('*.nev*');

%% load in cds and extract data

    cds = commonDataStructure();
    cds.file2cds([inputData.folderpath fileList(1).name],inputData.task,inputData.ranBy,...
        inputData.monkey,inputData.labnum,inputData.array1,inputData.mapFileName,'recoverPreSync');
    cd(pwd);
    
%% convert cds to trial data
    params.event_list = {'bumpTime';'bumpMagnitude';'bumpDir';'isStimTrial';'stimCode'};
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
    params.threshold_mult = 0.3;
    td_reward = getMoveOnsetAndPeak(td_reward,params);

%% plot a set of reaches aligned to go cue with reaction time markers
    opts.MAX_PLOT = 4;
    plotReachesTD(td_reward,opts);
    
%% plot reaction times for each cue
    opts = [];
    plotReactionTimeDataTD(td_reward,opts);

    
%% plot a psychometic curve fitting detection probability for the provided data
    
    
    
    
    