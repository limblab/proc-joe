%% determine filename and input data
    input_data.folderpath = 'C:\Users\jts3256\Desktop\Kramer_07312013_BD\';
%     input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\RetiredMonkeys\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
    
    input_data.task='taskBD';
    input_data.ranBy='ranByJoseph'; 
    input_data.array1='arrayLeftS1'; 
    input_data.monkey='monkeyKramer';
    input_data.labnum = 6;
    
    pwd=cd;
    cd(input_data.folderpath)
    fileList = dir('*s.nev*');
    cd(pwd)
    
    
%% load in cds and extract data
    td_all = [];

    for fileNumber = 1:numel(fileList)  
        cds = commonDataStructure();
        cds.file2cds([input_data.folderpath fileList(fileNumber).name],input_data.task,input_data.ranBy,...
            input_data.monkey,input_data.labnum,input_data.array1,input_data.mapFileName,'recoverPreSync');
        cd(pwd);


    end
    
    
%% convert to td
    params.event_list = {'tgtDir';'isPrimaryTgt';...
        'bumpTime';'bumpDir';'bumpMagnitude';'isStimTrial';'isPrimaryTgt'};
    params.trial_results = {'R','F'};
    params.extra_time = [1,2];
    params.exclude_units = [0,255];
    td_all = parseFileByTrial(cds,params);

    td_all = getMoveOnsetAndPeak(td_all);
    td_all = removeBadTrials(td_all);
    
%% split 

    params_split.split_idx_name = 'idx_bumpTime';
    params_split.extra_bins = [0,60];
    
    td_bump = splitTD(td_all,params_split);
    
    
%% prepare data for LDA
    meas = 

%% run LDA

    
    
    
    






