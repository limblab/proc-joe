%% set initial parameters

    input_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\Han_stim_evoked_move\Han_20200629_330Hz_chan16\';

%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\6250-002085\SN 6250-002085.cmp';
    
    input_data.date = '20190923';

    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyHan';
    input_data.ranBy = 'ranByJoe';
    input_data.lab = 6;
    input_data.mapFile = strcat('mapFile',mapFileName);
    input_data.task = 'taskCObump';

    pwd=cd;
    input_data.fileNum = 1;
  
    input_data.center_y = -33;
    input_data.center_x = 3;
    input_data.num_bins = 8;
    
%% make cds
    cd(input_data.folderpath)

    file_name = dir('*nev*');
    
    params.event_list = {'goCueTime';'bumpTime';'stimTime';'stimCode'};
    params.trial_results = {'R','A','I'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [255];

    move_onset_params.pre_move_thresh = 1000;
    move_onset_params.min_s = 3;
    move_onset_params.start_idx_offset = -10;
    move_onset_params.max_rt_offset = 400;
    td_all = [];
    for i = 1:numel(file_name)
        cds = commonDataStructure();
        cds.file2cds(strcat(input_data.folderpath,file_name(i).name),input_data.array,input_data.monkey,input_data.ranBy,...
            input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');

       
        td_temp = parseFileByTrial(cds,params);
        td_temp = stripSpikeSorting(td_temp);
        td_temp = getSpeed(td_temp);
    %     td_all = getMoveOnset(td_all,move_onset_params);
    %     td_all = removeBadTrials(td_all);
        td_all = [td_all,td_temp];
    end
    if(td_all(1).bin_size < 0.05)
        % set it to 50ms
        td_all = binTD(td_all,ceil(0.01/td_all(1).bin_size));
    end
    
        
%% get correlation during movement epoch
    corr_params = [];
    corr_params.signals = {'LeftS1_spikes'};
    td_move = trimTD(td_all,{'idx_goCueTime',0},{'idx_goCueTime',10});
    td_move = softNormalize(td_move);

%     avg_params = []; avg_params.conditions = 'all';
%     [td_move_avg,cond_idx] = trialAverage(td_move,{'target_direction'});
    td_move_mean_sub = subtractConditionMean(td_move);
    
    [rho,sort_idx] = pairwiseCorr(td_move_mean_sub,corr_params);


    
%% get PDs
    
    pd_params = [];
    pd_params.out_signals = 'LeftS1_spikes';
    pd_params.bootForTuning = 1;
    pd_params.num_boots = 100;
    pd_params.move_corr = 'vel';
    
    pd_all = getTDPDs(td_all,pd_params);
    
    
    
%% get PDs for bump

    td_bump = trimTD(td_all,{'idx_bumpTime',0},{'idx_bumpTime',300});

    pd_params = [];
    pd_params.out_signals = 'LeftS1_spikes';
    pd_params.bootForTuning = 0;
    pd_params.num_boots = 1;
    pd_params.move_corr = 'vel';
    
    pd_bump = getTDPDs(td_bump,pd_params);
    
%% get PDs for move
    td_move = trimTD(td_all,{'idx_goCueTime',0},{'idx_endTime',0});

    pd_params = [];
    pd_params.out_signals = 'LeftS1_spikes';
    pd_params.bootForTuning = 0;
    pd_params.num_boots = 1;
    pd_params.move_corr = 'vel';
    
    pd_move = getTDPDs(td_move,pd_params);

      
%% visualize
%% heatmap for preferred directions
disp('start')
    mapData = loadMapFile(mapFileName);
    
    optsPD.MAKE_BAR_PLOT = 1;

    optsPD.PLOT_CHANNELS = [1:96];
    optsPD.STIM_CHANNEL = 41;

    optsPD.MAX_RATIO = 1;
    optsPD.MIN_RATIO = -1;
    optsPD.LOG_SCALE = 0;
    optsPD.LOG_PARAM = 9;

    optsPD.FIGURE_SAVE = 0;
    optsPD.FIGURE_DIR = input_data.folderpath;
    optsPD.FIGURE_PREFIX = 'Han_20190924';
    
    [heatmapPD] = plotHeatmapsPD(td_all,pd_all,mapData,optsPD);
    
    
    
%% analyze movements during center hold caused by stim....
    window = [-20,50];
    td_stim = td_all(~isnan([td_all.idx_stimTime]) & [td_all.stimCode] == 2);
    td_stim = trimTD(td_stim,{'idx_stimTime',window(1)},{'idx_stimTime',window(2)});
    x_data = [window(1):1:window(2)]*td_stim(1).bin_size;
    % get correction latency....
    % find max speed within 600 ms from stim onset, then go backwards to
    % get movement onset
    move_on_idx = nan(numel(td_stim),1);
    center_hold_speed = nan(numel(td_stim),1);
    
    for i_trial = 1:numel(td_stim)
        [max_speed,max_speed_idx] = max(td_stim(i_trial).speed(td_stim(i_trial).idx_stimTime:end)); % relative to stim onset
        center_hold_speed(i_trial) = mean(td_stim(i_trial).speed(td_stim(i_trial).idx_stimTime-2:td_stim(i_trial).idx_stimTime+2));
        threshold = center_hold_speed(i_trial) + 0.1*(max_speed-center_hold_speed(i_trial));
        
        temp = find(td_stim(i_trial).speed(td_stim(i_trial).idx_stimTime:end) >= threshold,1,'first') - 1 + td_stim(i_trial).idx_stimTime - 1;
        if(~isempty(temp))
            move_on_idx(i_trial) = temp;
        end
    end
    
% make plots
    f_pos=figure(); hold on
    f_speed=figure(); hold on
    
    for i_trial = 1:numel(td_stim)
        % pos plot
        if(i_trial < 15)
            figure(f_pos);
            xlabel('X-pos (cm)');
            ylabel('Y-pos (cm)');
            plot(td_stim(i_trial).pos(1:td_stim(i_trial).idx_stimTime,1),td_stim(i_trial).pos(1:td_stim(i_trial).idx_stimTime,2),...
                'color',getColorFromList(1,1))
             plot(td_stim(i_trial).pos(td_stim(i_trial).idx_stimTime:end,1),td_stim(i_trial).pos(td_stim(i_trial).idx_stimTime:end,2),...
                'color',getColorFromList(1,0))
        end
        %speed plot
        if(center_hold_speed(i_trial) < 5)
            figure(f_speed)

            subplot(3,1,1); hold on
            ylabel('X-vel (cm/s)');
            plot(x_data,td_stim(i_trial).vel(:,1),'color',getColorFromList(1,1))
%             plot(x_data(move_on_idx(i_trial)),td_stim(i_trial).vel(move_on_idx(i_trial),1),'marker','.','markersize',20,'color',getColorFromList(1,2));
            
            subplot(3,1,2); hold on
            ylabel('Y-vel (cm/s)');
            plot(x_data,td_stim(i_trial).vel(:,2),'color',getColorFromList(1,1))
%             plot(x_data(move_on_idx(i_trial)),td_stim(i_trial).vel(move_on_idx(i_trial),2),'marker','.','markersize',20,'color',getColorFromList(1,2));

            subplot(3,1,3); hold on
            ylabel('Speed (cm/s)');
            xlabel('Time after stimulation onset (s)');
            plot(x_data,td_stim(i_trial).speed,'color',getColorFromList(1,1))
%             plot(x_data(move_on_idx(i_trial)),td_stim(i_trial).speed(move_on_idx(i_trial),1),'marker','.','markersize',20,'color',getColorFromList(1,2));
        end
    end


