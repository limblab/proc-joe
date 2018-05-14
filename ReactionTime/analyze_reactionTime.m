%% script to process reaction time data 
%% determine filename and input data
    inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Han_20180511_stim\';
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
td_all = [];
num_trials = 0;
for fileNumber = 1:numel(fileList)
    cds = commonDataStructure();
    cds.file2cds([inputData.folderpath fileList(fileNumber).name],inputData.task,inputData.ranBy,...
        inputData.monkey,inputData.labnum,inputData.array1,inputData.mapFileName,'recoverPreSync');
    cd(pwd);
    
    % convert cds to trial data
    params.event_list = {'isBumpTrial';'bumpTime';'bumpMagnitude';'bumpDir';'isStimTrial';'stimCode'};
    params.trial_results = {'R','F'};
    td_temp = parseFileByTrial(cds,params);
    td_temp = getGoCueTime(td_temp,cds);
    % append trial data into a single struct
    for t = 1:numel(td_temp)
        td_temp(t).trial_id = td_temp(t).trial_id + num_trials;
    end
    num_trials = num_trials + size(cds.trials,1);
    
    td_all = [td_all,td_temp];
end
%% separate out trials with results and go cue's
    [~, td_reward] = getTDidx(td_all, 'result', 'r');
    td_reward = td_reward(~isnan([td_reward.idx_goCueTime]));
    
%% get movement onset
    params.which_field = 'vel';
    params.field_idx = 1;
    params.start_idx_offset = 18;
    params.threshold_mult = 0.5;
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
    opts.WHICH_FIELD ='vel';
    opts.WHICH_IDX = [1];
    opts.BUMP_MAGS = [];
    opts.YLIM = [];
    plotReachesTD(td_reward,opts);
    
%% plot reaction times for each cue, psychometric curve
    td_reward_rt = td_reward(~isnan([td_reward.idx_movement_on]));
    opts = [];
    opts.FOLDER_PATH = inputData.folderpath;
    opts.SAVE_FIGURES = 0;
    opts.FIGURE_PREFIX = 'Han_20180509'; % no _ required
    opts.BUMP_MAGS = 0;
    opts.STIM_CODES = 0:7;
    opts.STIM_PARAMS = [50,100,150,200,250,300,350,400];

    
    opts.FIT = 1;
    opts.PLOT_BUMP = 1;
    opts.PLOT_STIM = 1;
    [data,plots] = plotReactionTimeDataTD(td_reward_rt,td_all,opts);


    
    
    
    