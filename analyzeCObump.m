%% set initial parameters

    input_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\CObump\Duncan_20190410\';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\6250-002085\SN 6250-002085.cmp';
    
    input_data.date = '20190410';
    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyDuncan';
    input_data.ranBy = 'ranByJoe';
    input_data.lab = 6;
    input_data.mapFile = strcat('mapFile',mapFileName);
    input_data.task = 'taskmulti_gadget';

    pwd=cd;
    input_data.fileNum = 1;
  
    input_data.center_y = -33;
    input_data.center_x = 3;
    input_data.num_bins = 8;
    
%%
    cd(input_data.folderpath)

    file_name = dir('*nev*');
    
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
    cds.file2cds(strcat(input_data.folderpath,file_name(1).name),input_data.array,input_data.monkey,input_data.ranBy,...
        input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
    
    removeIDFromUnits;
    
    %%
    td_all = parseFileByTrial(cds,params);
    td_all = getSpeed(td_all);
    td_all = removeBadTrials(td_all);
    td_all = getMoveOnset(td_all,move_onset_params);
    td_all = removeBadTrials(td_all);
    
%% get PDs
    
    pd_params = [];
    pd_params.out_signals = 'LeftS1_spikes';
    pd_params.bootForTuning = 0;
    pd_params.num_boots = 1;
    pd_params.move_corr = 'vel';
    
    pd_all = getTDPDs(td_all,pd_params);
    

    
    