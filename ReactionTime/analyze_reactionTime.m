%% script to process reaction time data 
%% determine filename and input data
    inputData.folderpath = 'C:\Users\jts3256\Desktop\Han_20180813\';
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
            inputData.monkey,inputData.labnum,inputData.array1,inputData.mapFileName);
        cd(pwd);

        if(~isempty(cds.trials) && size(cds.trials,1) > 1 && sum(cds.trials.result == 'R' | cds.trials.result == 'F') ~= 0)
            % convert cds to trial data
            params.event_list = {'bumpStaircaseIdx';'tgtOnTime';'isBumpTrial';'bumpTime';'bumpMagnitude';'bumpDir';'isStimTrial';'stimCode';'tgtDir'};
            params.trial_results = {'R','F'};
            params.extra_time = [0.8,2];
            params.include_ts = 1;
            params.exclude_units = [0,255];
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
%     clear cds
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
    opts.MAX_PLOT = 40;
    opts.WHICH_FIELD ='acc';
    opts.WHICH_IDX = [1];
    opts.BUMP_MAGS = [1];
    opts.YLIM = [];
    opts.STIM_CODES = [];
    plotReachesTD(td_reward,opts);
    
%% plot reaction times for each cue, psychometric curve
    opts = [];
    
    opts.SAVE_FIGURES = 0;

    td_reward_rt = td_reward(~isnan([td_reward.idx_movement_on]));
    opts.FOLDER_PATH = inputData.folderpath;
    opts.FIGURE_PREFIX = 'Han_20180809'; % no _ required
    opts.BUMP_MAGS = [];
    opts.STIM_CODES = [];
%     opts.STIM_PARAMS = [10,20,30,40,50,60,70,80,90,100];
%     opts.STIM_LABEL = 'Amplitude (\muA)';
%     opts.STIM_PARAMS = [50,100,150,200,250,300,350,400,450,500];
%     opts.STIM_LABEL = 'Frequency (Hz)';
%     opts.STIM_PARAMS = [25,50,75,100,125,150,175,200,250,300];
%     opts.STIM_LABEL = 'Train length (ms)';
    opts.STIM_PARAMS = [1,2,3,4,5,6,7,8,9];
    opts.STIM_LABEL = ['chan'];
    opts.FIT = 0;
    opts.PLOT_BUMP = 0;
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
    td_reward_rt = td_reward(~isnan([td_reward.idx_movement_on]));
    td_reward_rt_stim = td_reward_rt(isnan([td_reward_rt.idx_bumpTime]));
%     output_data = plotSpikesRT(td_reward_rt_stim,opts);
    output_data = cell(size(td_reward_rt_stim(1).LeftS1_spikes,2),1);
    rsquares = [];
    opts.TAU = 100000;
    for s = 1:size(td_reward_rt_stim(1).LeftS1_spikes,2)
        opts.STIM = 1;
        opts.SPIKE_LIST = s;
        opts.MAKE_PLOTS = 0;
        output_data{s} = plotSpikesRT(td_reward_rt_stim,opts);
%         stats = [output_data{s}.fits.stats];
        rsquares(s,:) = output_data{s}.corr_all_codes;
    end
    [~,sortIdx] = sort(rsquares,'ascend');
%%
corr_all = [];
for max_idx = 1
    opts.SPIKE_LIST = [];
    opts.TAU = 10000;
    opts.STIM = 1;
    opts.STIM_CODES = [];
    opts.MAKE_PLOTS = 1;
    corr_data = plotSpikesRT(td_reward_rt_stim,opts);
    corr_all(max_idx) = corr_data.corr_all_codes;
end
%% do a linear regression to predict reaction time from neural data
    td_reward_rt_stim = td_reward_rt(isnan([td_reward_rt.idx_bumpTime]));
    
    opts.SPIKE_LIST = [sortIdx(1:6)];
    opts.BIN_SIZE = 5;
    opts.WINDOW = [0,9];
    num_train = 100;
    train_temp = randperm(numel(td_reward_rt_stim));
    opts.TRAIN_IDX = train_temp(1:num_train);
    pred_data = predictRT(td_reward_rt_stim,opts);
%     pred_data.y_manual = -0.0035*(sum(pred_data.x(:,2:end),2)) + 0.32;
    test_idx = setdiff(1:numel(td_reward_rt_stim),opts.TRAIN_IDX);
    
    figure();
    plot(pred_data.y_true(opts.TRAIN_IDX),pred_data.y_pred(opts.TRAIN_IDX),'k.','markersize',12)
    hold on
    plot(0.01*[15:35],0.01*[15:35],'k--')
    plot(pred_data.y_true(test_idx),pred_data.y_pred(test_idx),'r.','markersize',12)
    sse = sum((pred_data.y_true - pred_data.y_pred).^2)
    
    
%% fitlm with stim code
p_vals = [];
mdl_rsquare = [];
for max_idx = 4
    opts.SPIKE_LIST = [sortIdx(1:max_idx)]; % 7 16 6 12 9 4 15 17
    opts.TAU = 0.2;
    opts.STIM_CODES = [];
    opts.MAKE_PLOTS = 0;
    corr_data = plotSpikesRT(td_reward_rt_stim,opts);

    rt = corr_data.rt; num_spikes = corr_data.num_spikes; code = corr_data.code;
    rt_adj = rt;
    code_unique = unique(code);
    mean_rt = zeros(numel(code_unique),1);
    for c = 1:numel(code_unique)
        mean_rt(c) = mean(rt(code == code_unique(c)));
        rt_adj(code == code_unique(c)) = rt_adj(code == code_unique(c)) - mean_rt(c);
    end
    
    tbl = table(rt,rt_adj,num_spikes,code);
    
    mdl_adj = fitlm(tbl,'rt_adj ~ num_spikes');
    mdl = fitlm(tbl,'rt ~ num_spikes + code');
    
    p_vals(end+1,:) = [mdl.Coefficients.pValue(2), mdl_adj.Coefficients.pValue(2)];
    mdl_rsquare(end+1,:) = [mdl.Rsquared.Ordinary, mdl_adj.Rsquared.Ordinary];
end

%% pick units based on distance from stimulated electrode
stim_elec = 20;
map_data = loadMapFile(inputData.mapFileName(8:end));
stim_elec_idx = find(map_data.chan == stim_elec);
stim_elec_pos = [map_data.row(stim_elec_idx), map_data.col(stim_elec_idx)];

unit_pos = [];
for unit_idx = 1:size(td_reward_rt_stim(1).LeftS1_unit_guide,1)
    map_data_idx = find(map_data.chan == td_reward_rt_stim(1).LeftS1_unit_guide(unit_idx,1));
    unit_pos(unit_idx,:) = [map_data.row(map_data_idx),map_data.col(map_data_idx)];
end

d = sqrt(sum((unit_pos-stim_elec_pos).^2,2));
[~,sortIdx] = sort(d,'ascend');
%% count # of spikes for each condition

stim_codes = unique(data_stim.code);
% bump_codes = unique(data_bump.code);

num_spikes_all = [];
labels = {};

for s = 1:numel(stim_codes)
    data_stim_mask = data_stim.code == stim_codes(s);
    num_spikes_all = [num_spikes_all; data_stim.num_spikes(data_stim_mask)];
    cueInfo_idx = find([data.cueInfo.stimCode] == stim_codes(s) & [data.cueInfo.bumpMag] == 0);
    data.cueInfo(cueInfo_idx).num_spikes = data_stim.num_spikes(data_stim_mask);
    for i = 1:sum(data_stim_mask)
        labels{end+1,1} = ['stim-',num2str(s)];
    end
end

% for b = 5:numel(bump_codes)
%     data_bump_mask = data_bump.code == bump_codes(b);
%     num_spikes_all = [num_spikes_all; data_bump.num_spikes(data_bump_mask)];
%     for i = 1:sum(data_bump_mask)
%         labels{end+1,1} = ['bump-',num2str(b)];
%     end
% end

figure();
boxplot(num_spikes_all,labels)

%%

figure();
for i = 1:7
    plot(mean(data.cueInfo(i).rt),mean(data.cueInfo(i).num_spikes),'.','color','k','markersize',12);
    hold on
end
    
       
    