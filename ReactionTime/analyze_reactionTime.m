%% script to process reaction time data 
%% determine filename and input data
    inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Han_20180510_stim\';
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
            params.event_list = {'bumpStaircaseIdx';'tgtOnTime';'isBumpTrial';'bumpTime';'bumpMagnitude';'bumpDir';'isStimTrial';'stimCode';'tgtDir'};
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
    opts.STIM_PARAMS = [10,20,30,40,50,60,70,80,90,100];
    opts.STIM_LABEL = 'Amplitude (\muA)';
%     opts.STIM_PARAMS = [50,100,150,200,250,300,350,400,450,500];
%     opts.STIM_LABEL = 'Frequency (Hz)';
%     opts.STIM_PARAMS = [25,50,75,100,125,150,200,250,300];
%     opts.STIM_LABEL = 'Train length (ms)';
    opts.FIT = 1;
    opts.PLOT_BUMP = 1;
    opts.PLOT_STIM = 1;
    [data,plots] = plotReactionTimeDataTD(td_reward_rt,td_all_rt,opts);

    data.opts = opts;

    
%% plot rasters
    opts = [];
    opts.X_LIMITS = [-0.3,0.5];
    opts.EXTRA_TIME = params.extra_time;
    opts.PLOT_RASTER = 1;
    opts.PLOT_PSTH = 0;
    opts.PLOT_BUMP = 0;
    opts.PLOT_STIM = 1;
    plotRasterPSTHRT(td_reward,opts);
    
%% plot # spikes vs rt for each trial
    td_reward_rt_stim = td_reward_rt(isnan([td_reward_rt.idx_bumpTime]));
    rt = zeros(numel(td_reward_rt_stim),1);
    num_spikes = zeros(numel(td_reward_rt_stim),1);
    offset = 12;
    tau = 0.2; % exponential decaying window
    exp_window = exp(-(mode([td_reward_rt_stim.bin_size])*(0:offset)')/tau);
    for t = 1:numel(td_reward_rt_stim)
        rt(t) = td_reward_rt_stim(t).bin_size*(td_reward_rt_stim(t).idx_movement_on - td_reward_rt_stim(t).idx_goCueTime);
        num_spikes(t) = sum(sum(td_reward_rt_stim(t).LeftS1_spikes(td_reward_rt_stim(t).idx_movement_on:td_reward_rt_stim(t).idx_movement_on+offset,:).*exp_window));
        stim_code(t) = td_reward_rt_stim(t).stimCode;
    end
    
    colors = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]/255;
    figure;
    for s = 1:numel(stim_codes)
        plot(rt(stim_code == stim_codes(s)),num_spikes(stim_code == stim_codes(s)),'.','markersize',12,'color',colors(s,:))
        hold on
    end
    l = legend('1','2','3','4','5','6','7','8','9');
    set(l,'box','off')
    [fitObj,stats] = fit(rt,num_spikes,'a*x+b');
%     for s = 1:numel(stim_codes)
%         figure();
%         plot(rt(stim_code == stim_codes(s)),num_spikes(stim_code == stim_codes(s)),'.','markersize',12,'color','k')
%         xlim([0.1,0.3]);
%     end
%% do a linear regression to predict reaction time from neural data? or kinematics?
    window = [0,10];
    for t = 1:numel(td_reward_rt)
        td_reward_rt(t).reaction_time = td_reward_rt(t).idx_movement_on - td_reward_rt(t).idx_goCueTime;
    end
    td_stim = td_reward_rt([td_reward_rt.isStimTrial]);
    td_stim = td_stim(~isnan(([td_stim.reaction_time])));
    td_stim = trimTD(td_stim,{'idx_goCueTime',window(1)},{'idx_goCueTime',window(2)});
    td_stim = binTD(td_stim,5);

    for t = 1:numel(td_stim)
        td_stim(t).LeftS1_spikes = reshape(td_stim(t).LeftS1_spikes,numel(td_stim(t).LeftS1_spikes),1);
    end
    mdl_struct = [];
    mdl_struct.model_type = 'linmodel';
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

    
    
    
    
    