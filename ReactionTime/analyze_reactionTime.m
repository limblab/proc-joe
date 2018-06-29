%% script to process reaction time data 
%% determine filename and input data
    inputData.folderpath = 'C:\Users\jts3256\Desktop\reactionTime_20180427\';
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
            params.event_list = {'bumpStaircaseIdx';'tgtOnTime';'isBumpTrial';'bumpTime';'bumpMagnitude';'bumpDir';'isStimTrial';'stimCode';'tgtDir'};
            params.trial_results = {'R','F'};
            params.extra_time = [1,2];
            params.include_ts = 0;
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
% sanitize td_all spike information since we are merging files, units go in and out
    % get master list of units
    master_list = [-1,-1];
    for trial = 1:length(td_all)
        for unit = 1:size(td_all(trial).LeftS1_unit_guide,1)
            % check if unit is in master_list
            master_idx = find(master_list(:,1)==td_all(trial).LeftS1_unit_guide(unit,1));
            if(isempty(master_idx) || (~isempty(master_idx) && sum(master_list(master_idx,2) == td_all(trial).LeftS1_unit_guide(unit,2))==0))
                master_list(end+1,:) = td_all(trial).LeftS1_unit_guide(unit,:);
            end
        end
    end
    master_list(1,:) = []; % remove dummy idx
    % adjust spike data to match master list of units
    for trial = 1:length(td_all)
        temp_spikes = zeros(size(td_all(trial).LeftS1_spikes,1),size(master_list,1));
        temp_ts = cell(size(master_list,1),1);
        for unit = 1:size(td_all(trial).LeftS1_unit_guide,1)
            master_idx = find(sum(master_list == td_all(trial).LeftS1_unit_guide(unit,:),2) == 2);
            temp_spikes(:,master_idx) = td_all(trial).LeftS1_spikes(:,unit);
            temp_ts{master_idx} = td_all(trial).LeftS1_ts{unit};
        end
        td_all(trial).LeftS1_unit_guide = master_list;
        td_all(trial).LeftS1_spikes = temp_spikes;
        td_all(trial).LeftS1_ts = temp_ts;
    end
%% separate out trials with results and go cue's
    [~, td_reward] = getTDidx(td_all, 'result', 'r');
    td_reward = td_reward(~isnan([td_reward.idx_goCueTime]));
    
%% get movement onset
    params.which_field = 'acc'; % peak acceleration
    params.field_idx = 1;
    params.start_idx_offset = 15;
%     params.threshold_mult = 0.4; % chosen so that bumps match with tactile reported in (Godlove, 2014). Probably shouldn't change
    params.be_aggressive = 1;
    params.threshold_mult = 2.5;
    params.min_s = 10;
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
    opts.WHICH_FIELD ='acc';
    opts.WHICH_IDX = [1];
    opts.BUMP_MAGS = [];
    opts.YLIM = [];
    opts.STIM_CODES = [4];
    plotReachesTD(td_reward,opts);
    
%% plot reaction times for each cue, psychometric curve
    opts = [];
    
    opts.SAVE_FIGURES = 0;

    td_reward_rt = td_reward(~isnan([td_reward.idx_movement_on]));
    opts.FOLDER_PATH = inputData.folderpath;
    opts.FIGURE_PREFIX = 'Han_20180528'; % no _ required
    opts.BUMP_MAGS = [];
    opts.STIM_CODES = [];
%     opts.STIM_PARAMS = [10,20,30,40,50,60,70,80,90,100];
%     opts.STIM_LABEL = 'Amplitude (\muA)';
%     opts.STIM_PARAMS = [50,100,150,200,250,300,350,400,450,500];
%     opts.STIM_LABEL = 'Frequency (Hz)';
    opts.STIM_PARAMS = [25,50,75,100,125,150,175,200,250,300];
    opts.STIM_LABEL = 'Train length (ms)';
    opts.FIT = 1;
    opts.PLOT_BUMP = 1;
    opts.PLOT_STIM = 1;
    [data,plots] = plotReactionTimeDataTD(td_reward_rt,td_all_rt,opts);

    data.opts = opts;

    
%% plot rasters
    opts = [];
    opts.X_LIMITS = [-0.3,0.5];
    opts.EXTRA_TIME = params.extra_time;
    opts.PLOT_RASTER = 0;
    opts.PLOT_PSTH = 1;
    opts.PLOT_BUMP = 0;
    opts.PLOT_STIM = 1;
    plotRasterPSTHRT(td_reward,opts);
    
%% plot # spikes vs rt for each trial
    td_reward_rt_stim = td_reward_rt(isnan([td_reward_rt.idx_bumpTime]));
%     output_data = plotSpikesRT(td_reward_rt_stim,opts);
    output_data = cell(size(td_reward_rt_stim(1).LeftS1_spikes,2),1);
    rsquares = [];
    opts.TAU = 0.2;
    for s = 1:size(td_reward_rt_stim(1).LeftS1_spikes,2)
        opts.SPIKE_LIST = s;
        output_data{s} = plotSpikesRT(td_reward_rt_stim,opts);
%         stats = [output_data{s}.fits.stats];
        rsquares(s,:) = output_data{s}.corr_all_codes;
    end
    
%%
    opts.SPIKE_LIST = []; % 15 16 20 13 1 8 5 7
    opts.TAU = 0.2;
    opts.STIM_CODES = [];
    corr_data = plotSpikesRT(td_reward_rt_stim,opts);
    corr_data.fit_all_codes.stats
    xlim([0.15,0.35]);
%     ylim([0,10]);
%% do a linear regression to predict reaction time from neural data
% fuck TD, write our own code
    td_reward_rt_stim = td_reward_rt(isnan([td_reward_rt.idx_bumpTime]));
    
    opts.SPIKE_LIST = [15 16 17 20];
    opts.BIN_SIZE = 5;
    opts.WINDOW = [0,14];
    num_train = 140;
    train_temp = randperm(numel(td_reward_rt_stim));
    opts.TRAIN_IDX = train_temp(1:num_train);
    pred_data = predictRT(td_reward_rt_stim,opts);
    test_idx = setdiff(1:numel(td_reward_rt_stim),opts.TRAIN_IDX);
    
    figure();
    plot(pred_data.y_true(opts.TRAIN_IDX),pred_data.y_pred(opts.TRAIN_IDX),'k.','markersize',12)
    hold on
    plot(0.01*[15:35],0.01*[15:35],'k--')
    plot(pred_data.y_true(test_idx),pred_data.y_pred(test_idx),'r.','markersize',12)


%% stepwise regression
    x = [];
    y = [];
    train_idx = randperm(numel(td_stim));
    train_idx = train_idx(1:floor(numel(train_idx)/2));
    test_idx = setdiff(1:numel(td_stim),train_idx);
    for t = train_idx
        x = [x;td_stim(t).LeftS1_spikes, td_stim(t).LeftS1_spikes_shift];
        y = [y;td_stim(t).reaction_time];
    end
    
    [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(x,y,'display','on','penter',0.2,'scale','on');
 % test model   
    b_use = b;
    b_use(~inmodel) = 0;
    intercept = stats.intercept;
    
    x_test = [];
    y_test = [];
    for t = test_idx
        x_test = [x;td_stim(t).LeftS1_spikes, td_stim(t).LeftS1_spikes_shift];
        y_test = [y;td_stim(t).reaction_time];
    end
    y_pred = intercept + x_test*b_use;
    figure();
    plot(y_test,y_pred,'.','markersize',12)
    hold on
    plot([15:30],[15:30],'k--')
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

    
        
%     window = [0,25];
%     mdl_struct = [];
% %     spike_list = [6,7,8,10,13,15,16,20];
%     spike_list = [15,20];
%     
%     for t = 1:numel(td_reward_rt)
%         td_reward_rt(t).reaction_time = td_reward_rt(t).idx_movement_on - td_reward_rt(t).idx_goCueTime;
%     end
%     td_stim = td_reward_rt([td_reward_rt.isStimTrial]);
%     td_stim = td_stim(~isnan(([td_stim.reaction_time])));
%     td_stim = trimTD(td_stim,{'idx_goCueTime',window(1)},{'idx_goCueTime',window(2)});
%     td_stim = binTD(td_stim,1);
%     num_shift = 0;
%     if(num_shift > 0)
%         p = {};
%         for i = 1:num_shift
%             p{end+1} = 'LeftS1_spikes';
%             p{end+1} = i;
%         end
%         td_stim = dupeAndShift(td_stim,p);
%         mdl_struct.in_signals = {'LeftS1_spikes',spike_list;'LeftS1_spikes_shift',spike_list};
%     else
%         mdl_struct.in_signals = {'LeftS1_spikes',spike_list};
%     end
%     td_stim = trimTD(td_stim,'start',1);
%     
%     mdl_struct.model_type = 'linmodel';
% 
%     mdl_struct.model_name = 'reaction_time_pred';
%     mdl_struct.glm_distribution = 'normal';
%     mdl_struct.out_signals = {'reaction_time','all'};
%     
%     mdl_struct.train_idx = 1:78;
%     
%     mdl_struct.do_lasso = false;
%     mdl_struct.lasso_lambda = 0.5;
%     mdl_struct.lasso_alpha = 1;
%     
%     [td_stim,mdl_info] = getModel(td_stim,mdl_struct);

    
    
    