%% script to process reaction time data 
%% determine filename and input data
    inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Han_20180521_stim\';
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

        if(~isempty(cds.trials) && size(cds.trials,1) > 1 && sum(cds.trials.result == 'R' | cds.trials.result == 'F') ~= 0)
            % convert cds to trial data
            params.event_list = {'tgtOnTime';'isBumpTrial';'bumpTime';'bumpMagnitude';'bumpDir';'isStimTrial';'stimCode';'tgtDir'};
            params.trial_results = {'R','F'};
            params.extra_time = [1,2];
            params.include_ts = 1;
            td_temp = parseFileByTrial(cds,params);
            td_temp = getGoCueTime(td_temp,cds);
            % append trial data into a single struct
            for t = 1:numel(td_temp)
                td_temp(t).trial_id = td_temp(t).trial_id + num_trials;
            end
            num_trials = num_trials + size(cds.trials,1);

            td_all = [td_all,td_temp];
        end
    end
    clear cds
    clear td_temp
%% separate out trials with results and go cue's
    [~, td_reward] = getTDidx(td_all, 'result', 'r');
    td_reward = td_reward(~isnan([td_reward.idx_goCueTime]));
    
%% get movement onset
    params.which_field = 'acc'; % peak acceleration
    params.field_idx = 1;
    params.start_idx_offset = 15;
    params.threshold_mult = 0.4; % chosen so that bumps match with tactile reported in (Godlove, 2014). Probably shouldn't change
    params.min_s = 50;
    params.pre_move_thresh = 50;
    params.be_aggressive = 1;
    td_reward = getMoveOnset(td_reward,params);
    
% put movement on back into td_all
    reward_idx = [td_reward.trial_id];
    td_all_rt = td_all;
    for td_reward_idx = 1:numel(td_reward)
        td_all_idx = find([td_all.trial_id] == td_reward(td_reward_idx).trial_id);
        td_all_rt(td_all_idx).idx_movement_on = td_reward(td_reward_idx).idx_movement_on;
    end
    clear td_all
%% plot a set of reaches aligned to go cue with reaction time markers
    opts.MAX_PLOT = 400;
    opts.WHICH_FIELD ='vel';
    opts.WHICH_IDX = [1];
    opts.BUMP_MAGS = [];
    opts.YLIM = [];
    opts.STIM_CODES = [];
    plotReachesTD(td_reward,opts);
    
%% plot reaction times for each cue, psychometric curve
    opts = [];
    
    opts.SAVE_FIGURES = 1;

    td_reward_rt = td_reward(~isnan([td_reward.idx_movement_on]));
    opts.FOLDER_PATH = inputData.folderpath;
    opts.FIGURE_PREFIX = 'Han_20180521'; % no _ required
    opts.BUMP_MAGS = [];
    opts.STIM_CODES = [];
    opts.STIM_PARAMS = [10,20,30,40,50,60,70,80,90,100];
    opts.STIM_LABEL = 'Amp (\muA)';
    
    opts.FIT = 1;
    opts.PLOT_BUMP = 1;
    opts.PLOT_STIM = 1;
    [data,plots] = plotReactionTimeDataTD(td_reward_rt,td_all_rt,opts);

    data.opts = opts;

    
%% plot rasters
    opts = [];
    opts.X_LIMITS = [-0.3,0.3];
    opts.EXTRA_TIME = params.extra_time;
    opts.PLOT_RASTER = 1;
    opts.PLOT_PSTH = 0;
    opts.PLOT_BUMP = 0;
    plotRasterPSTHRT(td_reward,opts);
    
    
%% do a glm to predict reaction time from neural data? or kinematics?
    window = [0,10];
    for t = 1:numel(td_reward)
        td_reward(t).reaction_time = td_reward(t).idx_movement_on - td_reward(t).idx_goCueTime;
    end
    td_stim = td_reward([td_reward.isStimTrial]);
    td_stim = td_stim(~isnan(([td_stim.reaction_time])));
    td_stim = binTD(td_stim,5);
    td_stim = trimTD(td_stim,{'idx_goCueTime'},{'idx_goCueTime',3});
    for t = 1:numel(td_stim)
        td_stim(t).LeftS1_spikes = td_stim(t).LeftS1_spikes(:)';
    end
    mdl_struct = [];
    mdl_struct.model_type = 'glm';
    mdl_struct.model_name = 'reaction_time_pred';
    mdl_struct.in_signals = {'LeftS1_spikes','all'};
    mdl_struct.out_signals = {'reaction_time','all'};
    [td_stim,mdl_info] = getModel(td_stim,mdl_struct);

%%
    
    
    
    
    % %% plot the reaction time to different days/parameter sets together based on
% % probability of detection
% 
%     % load in the two data files, feed them to the function below
%     
%     plotProbDetectRT({Han_20180515_data,Han_20180516_data},[]);
    
    % %% plot # spikes in a window vs. rt for each trial
%     rt = zeros(numel(td_reward),1)-1000;
%     num_spikes = zeros(numel(td_reward),1)-1000;
%     offset = [0,10]; % 10 ms bins, these are idxs
%     for trial = 1:numel(td_reward)
%         if(~isnan(td_reward(trial).idx_movement_on))
%             rt(trial) = (td_reward(trial).idx_movement_on - td_reward(trial).idx_goCueTime)*td_reward(trial).bin_size;
%             window = td_reward(trial).idx_goCueTime + offset;
%             num_spikes(trial) = sum(sum(td_reward(trial).LeftS1_spikes(window(1):window(2),:)));
%         end
%     end
%     
%     rt(rt < 0) = [];
%     num_spikes(num_spikes < 0) = [];
% %     num_spikes(rt > 0.35) = [];
% %     rt(rt > 0.35) = [];
%     
%     figure();
%     plot(rt,num_spikes,'k.','markersize',12)
%     [f,gof] = fit(rt,num_spikes,'a*x+b');
%     hold on
%     xData = min(rt):0.01:max(rt);
%     plot(xData,f.a*xData+f.b,'k--')
% %     xlim([0,0.5])

    
    
    
    
    