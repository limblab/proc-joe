%% script to process reaction time data 
%% determine filename and input data
    inputData.folderpath = 'D:\Lab\Data\ReactionTime\Han_20200929_FCreactTime_blocking\';
    
%     inputData.folderpath = 'D:\Lab\Data\ReactionTime\Han_20180427_training\';
    inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\right S1 20180919\SN 6251-001804.cmp';
    
    
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
            params.event_list = {'bumpStaircaseIdx';'tgtOnTime';'isBumpTrial';'bumpTime';'bumpMagnitude';'bumpDir';'isStimTrial';'stimCode';'tgtDir';'isVisualTrial'};
            params.trial_results = {'R','F'};
            params.extra_time = [1,2];
            params.include_ts = 0;

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

    clear td_temp
    if(numel(fileList) > 1)
    % sanitize td_all spike information since we are merging files, units go in and out
        % get master list of units
        master_list = [-1,-1];
        for trial = 1:length(td_all)
            for i_unit = 1:size(td_all(trial).LeftS1_unit_guide,1)
                % check if unit is in master_list
                master_idx = find(master_list(:,1)==td_all(trial).LeftS1_unit_guide(i_unit,1));
                if(isempty(master_idx) || (~isempty(master_idx) && sum(master_list(master_idx,2) == td_all(trial).LeftS1_unit_guide(i_unit,2))==0))
                    master_list(end+1,:) = td_all(trial).LeftS1_unit_guide(i_unit,:);
                end
            end
        end
        master_list(1,:) = []; % remove dummy idx
        % adjust spike data to match master list of units
        for trial = 1:length(td_all)
            temp_spikes = zeros(size(td_all(trial).LeftS1_spikes,1),size(master_list,1));
            temp_ts = cell(size(master_list,1),1);
            for i_unit = 1:size(td_all(trial).LeftS1_unit_guide,1)
                master_idx = find(sum(master_list == td_all(trial).LeftS1_unit_guide(i_unit,:),2) == 2);
                temp_spikes(:,master_idx) = td_all(trial).LeftS1_spikes(:,i_unit);
    %             temp_ts{master_idx} = td_all(trial).LeftS1_ts{unit};
            end
            td_all(trial).LeftS1_unit_guide = master_list;
            td_all(trial).LeftS1_spikes = temp_spikes;
    %         td_all(trial).LeftS1_ts = temp_ts;
        end
    end
    
       
%% separate out trials with results and go cue's
    [~, td_reward] = getTDidx(td_all, 'result', 'r');
    td_reward = td_reward(~isnan([td_reward.idx_goCueTime]));
    td_reward = getSpeed(td_reward);
    
    % get movement onset
    params.field_idx = 1;
    params.start_idx_offset = 100;
    params.be_aggressive = 1;
    params.which_field = 'speed';

    % Han's parameters
    params.threshold_acc = -1; % absolute threshold on acceleration, using this instead of threshold_mult
    params.min_s = 1;
    params.pre_move_thresh = 1;

%     Duncan's parameters
%     params.threshold_acc = 35;
%     params.pre_move_thresh = 50;
%     params.min_s = 100;
%     params.peak_idx_offset = [0,70];
%     params.max_rt_offset = 50;
%     
    
    params.use_emg = 0;
    params.emg_idx = 13;
    
    td_reward = getMoveOnset(td_reward,params);
    
    % put movement on back into td_all
    reward_idx = [td_reward.trial_id];
    td_all_rt = td_all;
    for td_reward_idx = 1:numel(td_reward)
        td_all_idx = find([td_all.trial_id] == td_reward(td_reward_idx).trial_id);
        td_all_rt(td_all_idx).idx_movement_on = td_reward(td_reward_idx).idx_movement_on;
    end
    
    
%% plot a set of reaches aligned to go cue with reaction time markers
    opts.MAX_PLOT = 30;
    opts.WHICH_FIELD = 'speed';
    opts.DIR = 90;
    
    opts.BUMP_MAGS = [];
    opts.STIM_CODES = [];
    opts.KEEP_ONLY_VISUAL_TRIALS = 0;
    opts.YLIM = [];
    opts.COLOR = 'k';%getColorFromList(1,1);
    opts.RANDOM = 1;
    
    plotReachesTD(td_reward,opts);
 
    
%% plot reaction times for each cue, psychometric curve
    opts = [];
    
    opts.SAVE_FIGURES = 0;
 
    td_reward_rt = td_reward(~isnan([td_reward.idx_movement_on]));

    opts.FOLDER_PATH = inputData.folderpath;

    opts.FIGURE_PREFIX = 'Han'; % no _ required
    
%     opts.BUMP_MAGS = [0.5:0.5:4.5];
    opts.STIM_CODES = [0,1];
    opts.STIM_LABEL = 'Channel';
    opts.STIM_PARAMS = [1,2];    
%     opts.STIM_X_LABEL = {'10','15','20','25','30','35'};
%     opts.STIM_PARAMS = [5:5:35];
%     opts.STIM_LABEL = 'Frequency (Hz)';
%     opts.STIM_PARAMS = [50:50:500];
%     opts.STIM_PARAMS = [25,50,75,100,125,150,200,250,300];
%     opts.STIM_LABEL = 'Train length (ms)';
%     opts.STIM_PARAMS = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8] + repmat([-0.1,0.1],1,8);
%     opts.STIM_X_LABEL = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'};
%     opts.STIM_COLOR_IDX = [0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7];
%     opts.STIM_COLOR_ALPHA = repmat([1,0.7],1,8);

%     opts.STIM_PARAMS = [4,4,4,6,6,6,8,8,8,12,12,12,24,24,24] + repmat([-0.25,0,0.25],1,5);
%     opts.STIM_X_LABEL = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'};
%     opts.STIM_COLOR_IDX = [1,0,2,1,0,2,1,0,2,1,0,2,1,0,2,1,0,2];
    
%     opts.STIM_PARAMS = [1:1:16];
%     opts.STIM_LABEL = 'Elec';
%     opts.COLOR_LIST = 2;
%     opts.STIM_LABEL = 'Bump Mag';
    opts.FIT = 	0;
    opts.PLOT_BUMP = 0;
    opts.PLOT_STIM = 1;
    
    opts.LINE_WIDTH = 1.5;
    [data,plots] = plotReactionTimeDataTD(td_reward_rt,td_all_rt,opts);

    data.opts = opts;

%% do stats
    % find max bumpMag id
    data.cueInfo(isnan([data.cueInfo.percent_respond])) = [];
    [~,max_bump_idx] = max([data.cueInfo.bumpMag]);
%     [~,max_stim_idx] = max([data.cueInfo.stimCode]);
%     max_bump_idx = 16;
    max_stim_idx = 4;
    rt_stim = data.cueInfo(1).rt;
    rt_bump = data.cueInfo(13).rt;
    
    tail = 'left';
%     if(mean(rt_stim) > mean(rt_bump))
%         tail = 'right';
%     else
%         tail = 'left';
%     end
    
%     [h,p,ci,stats] = ttest2(rt_stim,rt_bump,0.95,tail,'unequal');
    
    
[p,h,stats] = ranksum(rt_stim,rt_bump,'method','approximate','tail',tail)

    
    
    
%% analyze blocking experiment
    % get single and together (combined) trial datas
    offset = 0.025;
    single_cue_idx = [1:78,160:213];
    
    
    td_single = td_reward(single_cue_idx);
    td_combined = td_reward(setdiff([1:numel(td_reward)],single_cue_idx));
    
    td_single = td_single(~isnan([td_single.idx_movement_on]));
    td_combined = td_combined(~isnan([td_combined.idx_movement_on]));
    td_reward_rt = td_reward(~isnan([td_reward.idx_movement_on]));
    
    % get rt and code data
    rt_single = ([td_single.idx_movement_on] - [td_single.idx_goCueTime])*td_single(1).bin_size;
    rt_combined = ([td_combined.idx_movement_on] - [td_combined.idx_goCueTime])*td_combined(1).bin_size;
    rt_all = ([td_reward_rt.idx_movement_on] - [td_reward_rt.idx_goCueTime])*td_reward(1).bin_size;
    
    codes_single = [td_single.stimCode];
    codes_combined = [td_combined.stimCode];
    codes_all = [td_reward_rt.stimCode];
    
    unique_codes = unique(codes_single);
    
    % plot mean and std dev RT for each cue during single and combined
    % blocks
    figure(); hold on;
    for i_code = 1:numel(unique_codes)  
        % plot single
        scatter(unique_codes(i_code)-offset*ones(sum(codes_single==unique_codes(i_code)),1),rt_single(codes_single==unique_codes(i_code)),24,'markerfacecolor','k','markeredgecolor','none');
        alpha(0.3);
        x1=errorbar(unique_codes(i_code)-offset,mean(rt_single(codes_single==unique_codes(i_code))),std(rt_single(codes_single==unique_codes(i_code))),'k.','markersize',30,'linewidth',1);
        
        % plot combined
        scatter(unique_codes(i_code)+offset*ones(sum(codes_combined==unique_codes(i_code)),1),rt_combined(codes_combined==unique_codes(i_code)),24,'markerfacecolor','r','markeredgecolor','none');
        alpha(0.3);
        x2=errorbar(unique_codes(i_code)+offset,mean(rt_combined(codes_combined==unique_codes(i_code))),std(rt_combined(codes_combined==unique_codes(i_code))),'r.','markersize',30,'linewidth',1);
    end
    l=legend([x1,x2],'Alone','Together');
    set(l,'location','best');
    formatForLee(gcf);
    set(gca,'fontsize',14)
    xlabel('Code');
    ylabel('RT (s)');
    
    % plot RT for each cue throughout session, color each cue
    figure(); hold on;
    for i_code = 1:numel(unique_codes)
        trial_idx = find(codes_all == unique_codes(i_code));
        plot(trial_idx,rt_all(trial_idx),'.','color',getColorFromList(1,2+i_code),'markersize',20);
    end
    l=legend('Code 1','Code 2');
    set(l,'location','best');
    formatForLee(gcf);
    set(gca,'fontsize',14)
    xlabel('Trial');
    ylabel('RT (s)');
        
    
    % plot RT at beginning and end of block
    num_trials = 25;
    figure(); hold on;
    for i_code = 1:numel(unique_codes)
        % single
        code_mask = codes_single==unique_codes(i_code);
            % beginning of block
            beg_idx = find(code_mask,num_trials,'first');
            beg_mask = 1:1:beg_idx(end);
            rt = rt_single(code_mask(1:beg_idx(end)) & beg_mask);
            errorbar(i_code-1.5*offset,mean(rt),std(rt),'k.');
            
            % end of block
            end_idx = find(code_mask,num_trials,'last');
            end_mask = numel(code_mask):-1:end_idx(1);
            rt = rt_single(code_mask(end_idx(1):end) & end_mask);
            errorbar(i_code-.5*offset,mean(rt),std(rt),'k.');
        % combined
        code_mask = codes_combined==unique_codes(i_code);
            % beginning of block
            beg_idx = find(code_mask,num_trials,'first');
            beg_mask = 1:1:beg_idx(end);
            rt = rt_combined(code_mask(1:beg_idx(end)) & beg_mask);
            errorbar(i_code+.5*offset,mean(rt),std(rt),'k.');
            
            % end of block
            end_idx = find(code_mask,num_trials,'last');
            end_mask = numel(code_mask):-1:end_idx(1);
            rt = rt_combined(code_mask(end_idx(1):end) & end_mask);
            errorbar(i_code+1.5*offset,mean(rt),std(rt),'k.');
    end
        
    % for each cue during combined block, plot RT on repetitive trials vs.
    % alternating (preceding cue is the same or different)
    cue_repeats = [0,diff(codes_combined)==0]; % want second cue in sequence
    cue_changes = [0,diff(codes_combined)~=0];
     
    figure(); hold on
    for i_code = 1:numel(unique_codes)
        code_mask = codes_combined == unique_codes(i_code);
        % repeats
        x1=errorbar(unique_codes(i_code)-offset, mean(rt_combined(cue_repeats & code_mask)), std(rt_combined(cue_repeats & code_mask)),...
            'color',getColorFromList(1,0),'marker','.','markersize',20);
        % changes
        x2=errorbar(unique_codes(i_code)+offset, mean(rt_combined(cue_changes & code_mask)), std(rt_combined(cue_changes & code_mask)),...
            'color',getColorFromList(1,1),'marker','.','markersize',20);
    end
    l=legend([x1,x2],'repeats','changes');
    set(l,'location','best');
    xlabel('Code');
    ylabel('RT (s)');
    
%% get electrodes stimulated on each trial -- this is for the random elecs experiments
% need to feed in the data matrix, EL_all (list of electrodes),
% stim_code_all

% some matlab states have a variable named data that overrides the one I
% want. This causes an error and can be fixed by running code above again
    input_data = [];
    input_data.data = data;
    input_data.td_all = td_all;
    input_data.td_reward = td_reward_rt;
    input_data.EL_all = EL_all;
    input_data.stim_code_all = stim_code_all;
    input_data.map_file_name = inputData.mapFileName(8:end);
    
    opts = [];
    
    [electrode_list_data] = getElectrodesOnEachTrial(input_data,opts);
    
%% plot rt for each trial against fastest electrode in each list
% need to load in all_files_data
    % get min rt from electrode list and store
    electrode_list_data.min_rt = 1000+zeros(size(electrode_list_data.rt));

    for i = 1:numel(electrode_list_data.rt)
        for j = 1:numel(electrode_list_data.EL_list{i})
            idx = find(all_files_data.chan == electrode_list_data.EL_list{i}(j));
            if(all_files_data.mean_rt(idx) < electrode_list_data.min_rt(i))
                electrode_list_data.min_rt(i) = all_files_data.mean_rt(idx);
            end
        end
        if(electrode_list_data.min_rt == 1000)
            disp(num2str(i))
        end
    end
    
    f = figure();
    f.Name = 'Han_bestVsManyElectrodes';
    subplot(1,2,1)
    plot(electrode_list_data.rt,electrode_list_data.min_rt,'.','markersize',12,'color',getColorFromList(1,1));
    hold on
    xlim([0.1,0.3])
    ylim([0.1,0.3])
    plot([0,1],[0,1],'r--','linewidth',1.5)
    ylabel('Best individual electrode')
    xlabel('Many electrode')
    formatForLee(gcf)
    set(gca,'fontsize',14)

    subplot(1,2,2)
    bin_edges = -0.1:0.01:0.1;
    histogram(electrode_list_data.rt-electrode_list_data.min_rt,bin_edges);
    formatForLee(gcf)
    xlabel('Many - Best');
    ylabel('Number of trials');
    set(gca,'fontsize',14)

%     
%     opts.SAVE_FIGURES = 0;
%     opts.FOLDER_PATH = inputData.folderpath;
%     opts.FIGURE_PREFIX = 'Han_20181108'; % no _ required
    
    
%     [dist_data,plots] = plotElectrodeDistanceRT(input_data,opts);

    
    
    
%% use trial_id in td_all and td_reward, plus the entries in data to find the rt for each 
% electrode group







%% fix the code below so thats its easier to use...
    
      
%% plot rasters
% plot psth
    opts = [];
    opts.X_LIMITS = [-0.3,0.5];
    opts.EXTRA_TIME = params.extra_time;
    opts.PLOT_RASTER = 1;
    opts.PLOT_PSTH = 0;
    opts.PLOT_BUMP = 0;
    opts.PLOT_STIM = 1;
    plotRasterPSTHRT(td_reward,opts);
    
%% get number of spikes per trial and rt for each code
    td_reward_rt = td_reward(~isnan([td_reward.idx_movement_on]));
    td_reward_rt_stim = td_reward_rt(isnan([td_reward_rt.idx_bumpTime]));
    spike_window_t = [0,0.1]; % s, relative to stim onset
    baseline_window_t = [-0.4,-0.1]; % s, relative to stim onset
    
    % pre-allocate arrays, variables
    code_list = nan(size(td_reward_rt_stim));
    rt_list = nan(size(td_reward_rt_stim)); 
    chan_list = td_reward_rt_stim(1).LeftS1_unit_guide(:,1);
    num_spikes = nan(size(td_reward_rt_stim(1).LeftS1_spikes,2),size(td_reward_rt_stim,2)); % trials x unit
    baseline_spikes = nan(size(td_reward_rt_stim(1).LeftS1_spikes,2),1);
    
    spike_window_idx = floor(spike_window_t/td_reward_rt_stim(1).bin_size);
    baseline_window_idx = floor(baseline_window_t/td_reward_rt_stim(1).bin_size);
    correction_factor = diff(spike_window_t)/diff(baseline_window_t);
    
    for i_trial = 1:numel(td_reward_rt_stim)
        % get code
        code_list(i_trial) = td_reward_rt_stim(i_trial).stimCode;
        % get rt
        rt_list(i_trial) = (td_reward_rt_stim(i_trial).idx_movement_on - td_reward_rt_stim(i_trial).idx_goCueTime)*td_reward_rt_stim(i_trial).bin_size;
        % get num spikes in window
        idx_go_cue = td_reward_rt_stim(i_trial).idx_goCueTime;
        num_spikes(:,i_trial) = sum(td_reward_rt_stim(i_trial).LeftS1_spikes(idx_go_cue+spike_window_idx(1):idx_go_cue+spike_window_idx(2),:));
        % get baseline spike count (based on stim window size)
        baseline_spikes(:,i_trial) = correction_factor*sum(td_reward_rt_stim(i_trial).LeftS1_spikes(idx_go_cue+baseline_window_idx(1):idx_go_cue+baseline_window_idx(2),:));
    end
    
    spikes_above_baseline = num_spikes - baseline_spikes;
    stim_codes = unique(code_list);
    
%% plot # spikes per trial vs rt for each code
    
    for i_code = 1:numel(stim_codes)
        figure();
        code_mask = code_list == stim_codes(i_code);
        plot(rt_list(code_mask==1),sum(spikes_above_baseline(:,code_mask==1)),'.')
    end
    
    
    
%% plot mean # spikes vs mean rt for each code
    
    figure(); hold on;
    mean_rts = [];
    mean_spikes = [];
    
    for i_code = 1:numel(stim_codes)
        code_mask = code_list == stim_codes(i_code);
        plot(mean(rt_list(code_mask==1)),mean(sum(spikes_above_baseline(:,code_mask==1))),'k.','markersize',20)
        mean_rts(i_code) = mean(rt_list(code_mask==1));
        mean_spikes(i_code) = mean(sum(spikes_above_baseline(:,code_mask==1)));
    end
    
    
    
%% plot heatmap showing stim elec location and recorded FR
    map_data = loadMapFile(inputData.mapFileName(8:end));
%     stim_elec_list = [15,42,44,48,50];
    stim_elec_list = [62];
    
    for i_code=1:numel(stim_codes)
        figure();
        code_mask = code_list == stim_codes(i_code);
        val_plot = mean(spikes_above_baseline(:,code_mask==1)');
        val_plot(val_plot<0) = 0;
        val_plot = val_plot/max(val_plot);
        
        for i_unit = 1:numel(chan_list) % chan_list defined in the get FR portion
            map_idx = find(chan_list(i_unit) == map_data.chan);
            
            color_to_use = (val_plot(i_unit))*[1,0,0];
            rectangle('Position',[map_data.row(map_idx),map_data.col(map_idx),1,1],'FaceColor',color_to_use, 'EdgeColor','none');
            hold on
        end
        % plot stim chan
        map_idx = find(stim_elec_list==map_data.chan);
        rectangle('Position',[map_data.row(map_idx),map_data.col(map_idx),1,1],'FaceColor','none','EdgeColor','m');
        
        plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5)
        
        set(gca,'visible','off')
        xlim([1,11])
        ylim([1,11])
        axis square
    end
    



    
    
    
 %% get number of spikes per trial and result
    td_detect = td_all([td_all.result] == 'R' | [td_all.result] == 'F');
    td_detect = td_detect(~isnan([td_detect.idx_goCueTime]));
 
    spike_window_t = [0,0.1]; % s, relative to stim onset
    baseline_window_t = [-0.4,-0.1]; % s, relative to stim onset
    
    % pre-allocate arrays, variables
    code_list = nan(size(td_detect));
    is_detect_list = nan(size(td_detect)); 
    chan_list = td_detect(1).LeftS1_unit_guide(:,1);
    num_spikes = nan(size(td_detect(1).LeftS1_spikes,2),size(td_detect,2)); % trials x unit
    baseline_spikes = nan(size(td_detect(1).LeftS1_spikes,2),1);
    
    spike_window_idx = floor(spike_window_t/td_detect(1).bin_size);
    baseline_window_idx = floor(baseline_window_t/td_detect(1).bin_size);
    correction_factor = diff(spike_window_t)/diff(baseline_window_t);
    
    for i_trial = 1:numel(td_detect)
        % get code
        code_list(i_trial) = td_detect(i_trial).stimCode;
        % get rt
        is_detect_list(i_trial) = td_detect(i_trial).result == 'R';
        % get num spikes in windowtd_detect(i_trial)
        idx_go_cue = td_detect(i_trial).idx_goCueTime;
        num_spikes(:,i_trial) = sum(td_detect(i_trial).LeftS1_spikes(idx_go_cue+spike_window_idx(1):idx_go_cue+spike_window_idx(2),:));
        % get baseline spike count (based on stim window size)
        baseline_spikes(:,i_trial) = correction_factor*sum(td_detect(i_trial).LeftS1_spikes(idx_go_cue+baseline_window_idx(1):idx_go_cue+baseline_window_idx(2),:));
    end
    
    spikes_above_baseline = num_spikes - baseline_spikes;
    stim_codes = unique(code_list);
    stim_codes(isnan(stim_codes)) = [];
%% plot num spikes for each code with whether the trial was detected or not -- use neighborhood around stim elec
    % get distance mask first
    stim_elec = 62;
    neigh_dist = 500; % um
    
    
    distance_mask = ones(size(chan_list));
    map_data = loadMapFile(inputData.mapFileName(8:end));
    stim_chan_idx = find(map_data.chan == stim_elec);
    stim_chan_pos = [map_data.row(stim_chan_idx),map_data.col(stim_chan_idx)];
    
    for i_chan = 1:numel(chan_list)
        chan_idx = find(map_data.chan == chan_list(i_chan));
        chan_pos = [map_data.row(chan_idx),map_data.col(chan_idx)];
        distance_mask(i_chan) = 400*sqrt(sum((stim_chan_pos-chan_pos).^2)) < neigh_dist;
    end
    
    
    % plot sum of spikes for neurons within distance
    figure(); hold on;
    for i_code = 1:numel(stim_codes)
        code_mask = code_list == stim_codes(i_code);
        detected_mask = is_detect_list == 1;
        
        % plot undetected trials
        if(sum(code_mask & ~detected_mask) > 0)
            plot(i_code-0.1,sum(spikes_above_baseline(distance_mask==1,code_mask & ~detected_mask),1),'r.','markersize',10)
        end
        % plot detected trials
        if(sum(code_mask & detected_mask) > 0)
            plot(i_code+0.1,sum(spikes_above_baseline(distance_mask==1,code_mask & detected_mask),1),'b.','markersize',10)
        end
    end

 %% classifier to predict result
    code = 2;
    code_mask = code_list == code;
    
    results = is_detect_list(code_mask)';
    spikes = spikes_above_baseline(:,code_mask)';
        
    % K-fold cross val 
    cv = cvpartition(numel(results),'KFold',numel(results));
    
    mdl = fitcdiscr(spikes,results,'CVPartition',cv);
    pred = nan(cv.NumObservations,1);
    
    for i_set = 1:cv.NumTestSets
        % get testing dataset predictions
        test_mask = test(cv,i_set);
        
        pred(test_mask==1) = predict(mdl.Trained{i_set},spikes(test_mask==1,:));
    end

    percent_correct = sum(pred == results)/numel(results)
    sum(results==1)/numel(results)
    
    
 
    