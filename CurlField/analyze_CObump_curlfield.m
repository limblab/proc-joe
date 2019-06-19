%% set initial parameters

    input_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\CObump\Han_20160405_CObump\';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    
    input_data.date = '20160405';
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
    
%% load in td's
    cd(input_data.folderpath)

    baseline_file = dir('*baseline*nev*');
    adaptation_file = dir('*adaptation*nev*');
    washout_file = dir('*washout*nev*');
    
    params.event_list = {'goCueTime';'tgtDir'};
    params.trial_results = {'R'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [1:255];

    move_onset_params.pre_move_thresh = 1000;
    move_onset_params.min_s = 3;
    move_onset_params.start_idx_offset = -10;
    move_onset_params.max_rt_offset = 400;
    

    cds = commonDataStructure();
    cds.file2cds(strcat(input_data.folderpath,baseline_file.name),input_data.array,input_data.monkey,input_data.ranBy,...
        input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
    
    removeIDFromUnits;
    td_baseline = parseFileByTrial(cds,params);
    td_baseline = getSpeed(td_baseline);
    td_baseline = removeBadTrials(td_baseline);
    td_baseline = getMoveOnset(td_baseline,move_onset_params);
    td_baseline = removeBadTrials(td_baseline);
    
    
    cds = commonDataStructure();
    cds.file2cds(strcat(input_data.folderpath,adaptation_file.name),input_data.array,input_data.monkey,input_data.ranBy,...
        input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
    
    removeIDFromUnits;
    td_adapt = parseFileByTrial(cds,params);
    td_adapt = getSpeed(td_adapt);
    td_adapt = removeBadTrials(td_adapt);
    td_adapt = getMoveOnset(td_adapt,move_onset_params);
    td_adapt = removeBadTrials(td_adapt);
    
    cds = commonDataStructure();
    cds.file2cds(strcat(input_data.folderpath,washout_file.name),input_data.array,input_data.monkey,input_data.ranBy,...
        input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
    
    removeIDFromUnits;
    td_washout = parseFileByTrial(cds,params);
    td_washout = getSpeed(td_washout);
    td_washout = removeBadTrials(td_washout);
    td_washout = getMoveOnset(td_washout,move_onset_params);
    td_washout = removeBadTrials(td_washout);
    
    
    cd(pwd);
    
    
%% Plot reaches during each group (baseline, adapt, washout)
% we want example reaches from baseline, and then all reaches in adapt
% colored
    num_baseline_plot = 20;
    baseline_color = [0,0,0];
    adapt_color = [1,0,0];
    washout_color = [0,0,1];
    window = [0,0];
    
    figure();
    plotReaches(td_baseline,datasample([1:numel(td_baseline)],num_baseline_plot,'Replace',false),...
        'idx_goCueTime','idx_endTime',window,baseline_color)
    
    figure();
    plotReaches(td_adapt,[1:50],'idx_goCueTime','idx_endTime',window,adapt_color);
    
    figure();
    plotReaches(td_adapt,numel(td_adapt)+[-50:1:-1],'idx_goCueTime','idx_endTime',window,adapt_color);

    figure();
    plotReaches(td_washout,[1:50],'idx_goCueTime','idx_endTime',window,washout_color)
    
    
%% get and plot reach metrics (takeoff angle)
    baseline_color = [0,0,0];
    adapt_color = [1,0,0];
    washout_color = [0,0,1];
    
    base_metrics = getAdaptationMetrics(td_baseline);
        base_metrics.takeoff_angle = base_metrics.takeoff_angle - mean(base_metrics.takeoff_angle);
    adapt_metrics = getAdaptationMetrics(td_adapt);
        adapt_metrics.takeoff_angle = adapt_metrics.takeoff_angle - mean(base_metrics.takeoff_angle);
    washout_metrics = getAdaptationMetrics(td_washout);
        washout_metrics.takeoff_angle = washout_metrics.takeoff_angle - mean(base_metrics.takeoff_angle);
    
    figure();
    plot(1:numel(base_metrics.takeoff_angle),abs(base_metrics.takeoff_angle),'color',baseline_color);
    hold on
    plot(numel(base_metrics.takeoff_angle)-0.5+[0,0],[0,180],'--','color',[0.5,0.5,0.5],'linewidth',2)
    plot(numel(base_metrics.takeoff_angle)+(1:numel(td_adapt)),abs(adapt_metrics.takeoff_angle),'color',adapt_color);
    plot(numel(base_metrics.takeoff_angle)-0.5+numel(td_adapt)+[0,0],[0,180],'--','color',[0.5,0.5,0.5],'linewidth',2)
    plot(numel(base_metrics.takeoff_angle)+numel(adapt_metrics.takeoff_angle)+(1:numel(td_washout)),...
        abs(washout_metrics.takeoff_angle),'color',washout_color);
    
    ylim([0,45])
    
%% Tuning of single neurons
    pd_params = [];
    pd_params.out_signals = 'LeftS1_spikes';
    pd_params.bootForTuning = 1;
    pd_params.num_boots = 10;
    pd_params.move_corr = 'vel';
    
    pd_baseline = getTDPDs(td_baseline,pd_params);
    
    pd_params.trial_idx = numel(td_adapt) + [-numel(td_baseline):1:-1];
    pd_adapt = getTDPDs(td_adapt,pd_params);

%% plot pd_baseline vs. pd_adapt with conf bounds
    plot(pd_baseline.velPD,pd_adapt.velPD,'.','markersize',12)
    
    











    