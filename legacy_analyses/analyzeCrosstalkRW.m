 input_data.folderpath = 'D:\Lab\Data\DLC_videos\Rocket_20210709_area2_rw\'; % DLC project folder
    
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Crackle 18E2\Map Files\Left S1\SN 6251-001695.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Rocket_19L1\MapFile\3a-area2\6251-002088\SN 6251-002088.cmp';
    
    use_td = 1;
    
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

    nev_file_name = dir('*nev*');
        
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
            curr_task = 'tasknone';
        else
            curr_task = input_data.task;
        end
        cds = commonDataStructure();
        cds.file2cds(strcat(nev_file_name(i_file).folder,filesep,nev_file_name(i_file).name),input_data.array,input_data.monkey,input_data.ranBy,...
            input_data.lab,input_data.mapFile,curr_task,'recoverPreSync','unsanitizedTimes');
        
        cds_list{i_file} = cds;
        clear cds;
        
    end
    % sort list based on order experiment was done (num_list)
    [num_list, sort_idx] = sort(num_list);
    task_list = task_list(sort_idx);
    cds_list = cds_list(sort_idx);
    nev_file_name = nev_file_name(sort_idx);
    
%% make trial data, strip spike sorting

    td_params = [];
    td_params.include_ts = 0;
    td_params.exclude_units = [255];
    
    td_params.noTrials=false;
    for i_cds = 1:numel(cds_list)
        td_list{i_cds} = parseFileByTrial(cds_list{i_cds},td_params);
        if(isempty(strfind(task_list{i_cds},'RT3D')) || isempty(cds_list{i_cds}.trials))
            td_list{i_cds} = removeTrials(td_list{i_cds});
            td_list{i_cds} = binTD(td_list{i_cds},round(0.01/td_list{i_cds}(1).bin_size));
        end
        td_list{i_cds} = stripSpikeSorting(td_list{i_cds});
    end

    
%% cross corr across the array
    
    tree_perm_use = []; max_corr = -100;
        % get td from td_list based on task_list
    space_idx = 1;
    td_use = td_list{space_idx};
        
    % get correlations across neurons
    corr_mat = corr(td_use.LeftS1_spikes);
        
    % cluster correlations  if correct task
    T = linkage(corr_mat);
    f_to_close = figure();
    [~,~,tree_perm_use] = dendrogram(T,100);
    close(f_to_close);
    % get max corr for color limits later
    for i = 1:size(corr_mat,1) corr_mat(i,i) = nan; end % set diagonal as nan
    max_corr = max(max_corr,max(max(corr_mat)));
        
    % plot correlations
    f=figure();
    corr_mat_use = corr(td_use.LeftS1_spikes(:,tree_perm_use));
    for i = 1:size(corr_mat_use,1) corr_mat_use(i,i) = nan; end % set diagonal as nan
    imagesc(corr_mat_use);
    b=colorbar;
    caxis(max_corr*[-1,1]);
    cmap = colorcet('D1'); colormap(cmap);

%% get max cross correlation for each neuron, plot as array map

    corr_map = zeros(10,10);
    corr_mat = corr(td_use.LeftS1_spikes);
    for i = 1:size(corr_mat,1) corr_mat(i,i) = -1000; end % set diagonal as nan
    map_data = loadMapFile(mapFileName);
    
    for i=1:size(corr_mat_use,1)
        max_corr = max(corr_mat(i,:));
        
        % get row and col for this elec
        map_idx = find(td_use.LeftS1_unit_guide(i,1) == map_data.chan);
        pos = [11-map_data.row(map_idx), map_data.col(map_idx)];
        
        corr_map(pos(1),pos(2)) = max_corr;
    end
    
    figure();
    imagesc(corr_map);
    b=colorbar;
    caxis([0,1]);
    cmap = colorcet('L1'); colormap(cmap);



    