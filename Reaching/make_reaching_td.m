%% setup initial parameters and data folder   

    input_data.folderpath = 'D:\Lab\Data\DLC_videos\Han_20201204_rwFreeReach\'; % DLC project folder
    
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    
    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyHan';
    input_data.ranBy = 'ranByJoe';
    input_data.lab = 6;
    input_data.mapFile = strcat('mapFile',mapFileName);
    input_data.task = 'taskRW';

    pwd=cd;
    input_data.fileNum = 1;
  
    input_data.center_y = -33;
    input_data.center_x = 3;
    input_data.num_bins = 8;
    
    
%% make cds
    cd(input_data.folderpath)

    nev_file_name = dir('neural-data\*nev*');
    
    % remove filename with .mat, also switch .mat to .nev
    keep_mask = ones(size(nev_file_name));
    for i_nev = 1:numel(nev_file_name)
        nev_file_name(i_nev).name = [nev_file_name(i_nev).name(1:end-3),'mat'];
        if(~isempty(strfind(nev_file_name(i_nev).name,'-s')))
            keep_mask(i_nev) = 0;
        end
    end
    nev_file_name = nev_file_name(keep_mask==1);
    dlc_file_name = dir('reconstructed-3d-data\*.csv');
    
    num_list = nan(numel(nev_file_name),1);
    task_list = cell(numel(nev_file_name),1);
    cds_list = cell(numel(nev_file_name),1);
    td_list = cell(numel(nev_file_name),1);
    
    for i_file = 1:numel(nev_file_name)
        % number is last 3 values in filename
        num_list(i_file) = str2num(nev_file_name(i_file).name(end-6:end-4));
                
        % task is between 2nd and 3rd underscore
        underscore_idx = strfind(nev_file_name(i_file).name,'_');
        task_list{i_file} = nev_file_name(i_file).name(underscore_idx(2)+1:underscore_idx(3)-1);
        
        % load in cds, save extension so that we can line 3D tracking data
        % up with correct file. Also save which experiment type this is
        if(~isempty(strfind(task_list{i_file},'RT3D')))
            curr_task = 'taskUnknown';
        else
            curr_task = input_data.task;
        end
        cds = commonDataStructure();
        cds.file2cds(strcat(nev_file_name(i_file).folder,filesep,nev_file_name(i_file).name),input_data.array,input_data.monkey,input_data.ranBy,...
            input_data.lab,input_data.mapFile,curr_task,'recoverPreSync','ignoreJumps','ignoreFilecat');
        
        cds_list{i_file} = cds;
        clear cds;
        
    end
    % sort list based on order experiment was done (num_list)
    [num_list, sort_idx] = sort(num_list);
    task_list = task_list(sort_idx);
    cds_list = cds_list(sort_idx);
    nev_file_name = nev_file_name(sort_idx);
    
    
% load in 3D reaching data and place in cds. Need entire num_list since that is the order the csv files are in 
    % e.g: csv_0 = min value in num_list, then it increases and so on
    for i_file = 1:numel(dlc_file_name)
        cds_list{i_file}.loadRawMarkerDataDLC([dlc_file_name(i_file).folder,filesep,dlc_file_name(i_file).name]);
    end

%% convert to trial data
    td_params = [];
    td_params.all_points = 1;
    td_params.include_ts = 1;
    for i_cds = 1%:numel(cds_list)
        if(~isempty(strfind(task_list{i_cds},'RT3D')))
            td_params.noTrials=true;
        else
            td_params.noTrials=false;
        end
        td_list{i_cds} = parseFileByTrial(cds_list{i_cds},td_params);
        if(isempty(strfind(task_list{i_cds},'RT3D')))
            td_list{i_cds} = removeTrials(td_list{i_cds});
            td_list{i_cds} = binTD(td_list{i_cds},round(0.01/td_list{i_cds}(1).bin_size));
        end
        
    end

