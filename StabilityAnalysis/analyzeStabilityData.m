%% set initial parameters

    input_data.folderpath = 'C:\Users\jts3256\Desktop\Duncan_stim_data\Duncan_CObump\';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\L_MillerLab\limblab-archive\Retired Animal Logs\Monkeys\Kramer 10I1\Kramer sept 2012 implant array mapping\6251-0922.cmp';
%     mapFileName = '\\fsmresfiles.fsm.northwestern.edu\fsmresfiles\Basic_Sciences\Phys\L_MillerLab\limblab-archive\Retired Animal Logs\Monkeys\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
    
    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyDuncan';
    input_data.ranBy = 'ranByJoe';
    input_data.lab = 6;
    input_data.mapFile = strcat('mapFile',mapFileName);
    input_data.task = 'taskCObump';

    pwd=cd;
  
    input_data.center_y = -33;
    input_data.center_x = 3;
    input_data.num_bins = 8;
    
%% make trial data for each file
    cd(input_data.folderpath)
    folders = dir(['*',input_data.monkey(7:end),'*']);
    
    params.event_list = {'goCueTime';'bumpTime'};
    params.trial_results = {'R'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [255];

    move_onset_params.pre_move_thresh = 1000;
    move_onset_params.min_s = 3;
    move_onset_params.start_idx_offset = -10;
    move_onset_params.max_rt_offset = 400;
    
    td_all_files = cell(numel(folders),1);
    
    for folder_idx = 1%:numel(folders)
        td_all = [];
        cd([folders(folder_idx).folder,'\',folders(folder_idx).name]);
        file_name = dir('*nev*');

        for file_idx = 1:numel(file_name)
            disp(file_name(file_idx).name);
                cds = commonDataStructure();
                cds.file2cds(strcat(file_name(file_idx).folder,'/',file_name(file_idx).name),input_data.array,input_data.monkey,input_data.ranBy,...
                    input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');

    %             % REMOVE ID FROM UNITS and make trial data
                td_temp = parseFileByTrial(cds,params);
                td_temp = stripSpikeSorting(td_temp);
                td_temp = getSpeed(td_temp);
                td_temp = binTD(td_temp,50);
                td_all = [td_all, td_temp];
        end
        
        td_all_files{folder_idx} = td_all;

    end

%% remove move bump trials
    if(strcmpi(input_data.task,'taskCObump'))
        td_all_use = td_all_files;
        for t = 1:numel(td_all_files)
            keep_mask = isnan([td_all_use{t}.idx_bumpTime]);
            td_all_use{t} = td_all_use{t}(keep_mask);
        end
    else
        td_all_use = td_all_files;
    end
    keep_mask = ~cellfun(@isempty,td_all_use);
    td_all_use = td_all_use(keep_mask);

%% get PDs across whole file
% get normalized firing rate during a reach in each direction for each file
    pd_all = zeros(130,numel(td_all_use)); % more rows than required
    mod_depth_all = zeros(size(pd_all));
    pd_all_ci = zeros(130,numel(td_all_use),2);
    mod_depth_all_ci = zeros(size(pd_all_ci));
    is_vel_tuned = zeros(size(pd_all));
    
    for t = 1:numel(td_all_use)
        pd_params = [];
        pd_params.out_signals = 'LeftS1_spikes';
        pd_params.bootForTuning = 1;
        pd_params.num_boots = 100;
            pd_params.move_corr = 'vel';

        pd_temp = getTDPDs(td_all_use{t},pd_params);
        pd_all(td_all_use{t}(1).LeftS1_unit_guide(:,1),t) = pd_temp.velPD;
        mod_depth_all(td_all_use{t}(1).LeftS1_unit_guide(:,1),t) = pd_temp.velModdepth;
        pd_all_ci(td_all_use{t}(1).LeftS1_unit_guide(:,1),t,:) = pd_temp.velPDCI;
        mod_depth_all_ci(td_all_use{t}(1).LeftS1_unit_guide(:,1),t,:) = pd_temp.velModdepthCI;
        is_vel_tuned(td_all_use{t}(1).LeftS1_unit_guide(:,1),t) = pd_temp.velTuned;
        
    end
  
    if(strcmpi(input_data.task,'CObump'))
        target_dirs = unique([td_all_use{1}.target_direction]);
        target_dirs = target_dirs(~isnan(target_dirs));
        z_score_all = cell(numel(target_dirs),1); % each index is a reach direction
    end
    
%% get z_score based on either bumps or reaches in the same direction
    if(strcmpi(input_data.task,'CObump'))
        for t = 1:numel(td_all_use)
            for tgt_dir_idx = 1:numel(target_dirs)
                td_tgt_dir = td_all_use{t}([td_all_use{t}.target_direction] == target_dirs(tgt_dir_idx));

                td_tgt_dir = binTD(td_tgt_dir,50); % go to 50 ms bin
                td_tgt_dir = td_tgt_dir(~isnan([td_tgt_dir.idx_goCueTime]) & ~isnan([td_tgt_dir.idx_endTime]));

                td_tgt_dir_base = trimTD(td_tgt_dir,{'idx_goCueTime',-10},{'idx_goCueTime',-2});
                td_tgt_dir_base = joinTrials(td_tgt_dir_base);

                td_tgt_dir_move = trimTD(td_tgt_dir,{'idx_goCueTime',0},{'idx_endTime',0});
                td_tgt_dir_move = joinTrials(td_tgt_dir_move);

                mean_fr_move = mean(td_tgt_dir_move.LeftS1_spikes)/td_tgt_dir_move.bin_size;
                std_fr_base = std(td_tgt_dir_base.LeftS1_spikes)/td_tgt_dir_base.bin_size;
                mean_fr_base = mean(td_tgt_dir_base.LeftS1_spikes)/td_tgt_dir_base.bin_size;

                z_score_dir = ((mean_fr_move - mean_fr_base)./std_fr_base)'; % make column vector
                z_score_all{tgt_dir_idx}(:,t) = z_score_dir;
            end
        end
    end

    
%% get change in PD, mod depth, and z_score for each direction across files and visualize
    
    file_idx = figure();
    chan_list = [30,63,68,90,52,72,73,96,3,25,27,71,26,33,66,70,29,47,52,75,6,13,23,62,3,29,88,95,46,52,55,65,18,24,67,70];
    plot(pd_all(chan_list,:)'*180/pi - pd_all(chan_list,1)'*180/pi);

%     f = figure();
% 
%     plot(mod_depth_all');
%     
%     for i = 1:numel(z_score_all)
%         f = figure();
%         plot(z_score_all{i}');
%     end

%% histogram of differences
    compare_idx = [1,3];
%     chan_list = [13,15,27,61,24,26,32,95,4,5,29,91,1,2,9,10,7,11,44,64,9,19,62,4,25,85,96,8,23,44,60,22,62];
%     chan_list = [3,43,49,70,35,44,66,87,46,48,50,51,2,33,67,74];
    keep_mask = is_vel_tuned(1:96,compare_idx(1)) & is_vel_tuned(1:96,compare_idx(2));
%     keep_mask = zeros(96,1); keep_mask(chan_list) = 1;
  
    figure();
    histogram(angleDiff(pd_all(keep_mask==1,compare_idx(1)),pd_all(keep_mask==1,compare_idx(2)),1),20)


    
    
    
