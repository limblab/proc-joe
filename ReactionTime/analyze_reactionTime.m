%% script to process reaction time data 
%% determine filename and input data
    inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Han_20180509_stimTraining\';
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
    params.which_field = 'acc'; % peak acceleration
    params.field_idx = 1;
    params.start_idx_offset = 15;
    params.threshold_mult = 0.4; % chosen so that bumps match with tactile reported in (Godlove, 2014). Probably shouldn't change
%     params.min_s = 10;
    params.min_s = 50;
    params.pre_move_thresh = 50;
    td_reward = getMoveOnset(td_reward,params);
    
% put movement on back into td_all
    reward_idx = [td_reward.trial_id];
    td_all_rt = td_all;
    for td_reward_idx = 1:numel(td_reward)
        td_all_idx = find([td_all.trial_id] == td_reward(td_reward_idx).trial_id);
        td_all_rt(td_all_idx).idx_movement_on = td_reward(td_reward_idx).idx_movement_on;
    end
            
%% plot a set of reaches aligned to go cue with reaction time markers
    opts.MAX_PLOT = 20;
    opts.WHICH_FIELD ='vel';
    opts.WHICH_IDX = [1];
    opts.BUMP_MAGS = [];
    opts.YLIM = [];
    opts.STIM_CODES = [0:6]';
    plotReachesTD(td_reward,opts);
    
%% plot reaction times for each cue, psychometric curve
    opts = [];
    
    opts.SAVE_FIGURES = 1;

    td_reward_rt = td_reward(~isnan([td_reward.idx_movement_on]));
    opts.FOLDER_PATH = inputData.folderpath;
    opts.FIGURE_PREFIX = 'Han_20180509'; % no _ required
    opts.BUMP_MAGS = [];
    opts.STIM_CODES = [];
    opts.STIM_PARAMS = [40,50,60,70,80,90,100];
    opts.STIM_LABEL = 'Amplitude (\muA)';
    
    opts.FIT = 1;
    opts.PLOT_BUMP = 1;
    opts.PLOT_STIM = 1;
    [data,plots] = plotReactionTimeDataTD(td_reward_rt,td_all_rt,opts);

    data.opts = opts;
%% plot the reaction time to different days/parameter sets together based on
% probability of detection

    % load in the two data files, feed them to the function below
    
    plotProbDetectRT({Han_20180510_data,Han_20180511_data},[]);
    
    