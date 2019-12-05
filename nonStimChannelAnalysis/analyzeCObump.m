%% set initial parameters

    input_data.folderpath = 'E:\Data\Joseph\Duncan_stim_data\Duncan_20191025_CObump\';

    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'Z:\Basic_Sciences\Phys\L_MillerLab\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\6250-002085\SN 6250-002085.cmp';
    
    input_data.date = '20191025';
    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyDuncan';
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
    
    params.event_list = {'goCueTime';'bumpTime'};
    params.trial_results = {'R'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [1:255];

    move_onset_params.pre_move_thresh = 1000;
    move_onset_params.min_s = 3;
    move_onset_params.start_idx_offset = -10;
    move_onset_params.max_rt_offset = 400;
    

    cds = commonDataStructure();
    cds.file2cds(strcat(input_data.folderpath,file_name(1).name),input_data.array,input_data.monkey,input_data.ranBy,...
        input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
    
    
%% REMOVE ID FROM UNITS and make trial data
   
    td_all = parseFileByTrial(cds,params);
    td_all = stripSpikeSorting(td_all);
    td_all = getSpeed(td_all);
    td_all = removeBadTrials(td_all);
%     td_all = getMoveOnset(td_all,move_onset_params);
%     td_all = removeBadTrials(td_all);

    if(td_all(1).bin_size < 0.01)
        % set it to 50ms
        td_all = binTD(td_all,ceil(0.05/td_all(1).bin_size));
    end
    
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
%% get PD data in heatmap form
    mapData = loadMapFile(mapFileName);
    [pdHeatmapData,pdAlphaData] = getHeatmapDataPD(td_all,pd_all,mapData); 

%% plot heatmaps
    
    optsPD.MAKE_BAR_PLOT = 1;
    
%     plot_chans = [1:96];
%     plot_chans = cell2mat(arrayData{1}.CHAN_LIST);
    plot_chans = [];
    for i = 1:numel(arrayData)
        plot_chans(i) = arrayData{i}.CHAN_REC;
    end
    optsPD.PLOT_CHANNELS = plot_chans;%[1:96];
    
    optsPD.FIGURE_SAVE = 0;
    optsPD.FIGURE_DIR = input_data.folderpath;
    optsPD.FIGURE_PREFIX = 'Han_20190924';
    
    [heatmapPD,pdHeatmapData,alphaData] = plotHeatmapsPD(td_all,pd_all,mapData,optsPD);
    
    figure
    histogram(pdHeatmapData(alphaData == 1),18)
    [PDBinCount,PDBinEdges] = histcounts(pdHeatmapData(alphaData == 1),18);
    