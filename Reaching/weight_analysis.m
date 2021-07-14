%% Set up meta info and load trial data
    clear; clc;
    if ispc
        folderpath = 'D:\Lab\Data\DLC_videos\Han_20210709_freeReachWeight\neural-data\';
    else
        folderpath = '/data/raeed/project-data/limblab/s1-kinematics';
    end
    
    % load data
    file_info = dir(fullfile(folderpath,'*td*'));

    filenames = horzcat({file_info.name})';
    
    % save directory information (for convenience, since this code takes a while)
    savefile = true;
    if savefile
        savedir = fullfile(folderpath,'reaching_experiments','EncodingResults');
        if ~exist(savedir,'dir')
            mkdir(savedir)
        end
        run_date = char(datetime('today','format','yyyyMMdd'));
        savename = sprintf('encoderResults_run%s.mat',run_date);
    end
    
    arrayname = 'LeftS1';
    bin_size = 0.05; % s
    

%% Loop through trial data files to clean up 
    td_all = {};
    task_list_all = {};
    robot_height_all = {};
    for filenum = 1:length(filenames)
        % Load data
        load(fullfile(file_info(filenum).folder,file_info(filenum).name));
        keep_mask = ones(size(td_list));
        robot_height = [];
        for i_td = 1:numel(td_list) % each entry in td_list is a different trial_data for a different experiment (free reach vs 2D random walk for example)
            % resample trial data to appropriate bin size
            if(td_list{i_td}.bin_size <= bin_size)
                td_list{i_td} = binTD(td_list{i_td},bin_size/td_list{i_td}.bin_size);
            else
                warning('td bin size is larger than desired bin size');
            end

            % remove unsorted neurons
            unit_ids = td_list{i_td}(1).(sprintf('%s_unit_guide',arrayname));
            unsorted_units = (unit_ids(:,2)==0);
            new_unit_guide = unit_ids(~unsorted_units,:);
            for trialnum = 1:length(td_list{i_td})
                td_list{i_td}(trialnum).(sprintf('%s_unit_guide',arrayname)) = new_unit_guide;

                spikes = td_list{i_td}(trialnum).(sprintf('%s_spikes',arrayname));
                spikes(:,unsorted_units) = [];
                td_list{i_td}(trialnum).(sprintf('%s_spikes',arrayname)) = spikes;
            end
            
            % add firing rates in addition to spike counts
            td_list{i_td} = addFiringRates(td_list{i_td},struct('array',arrayname));

            if(isfield(td_list{i_td},'dlc_pos'))
                % set origin as shoulder position at t=0
                td_list{i_td} = setOriginAsShoulder(td_list{i_td},0); % use fixed position (t=0) or set shoulder as 0 for each data point.
                % get marker velocity
                td_list{i_td} = getDifferential(td_list{i_td},struct('signals','dlc_pos','alias','dlc_vel'));                
            end
            % get trial start and end based on button.
            [trial_start,trial_end, experiment_phase_mask] = getButtonHigh(td_list{i_td},'ana_var');
            
        end
        td_all{filenum} = td_list;
        task_list_all{filenum} = task_list;
        robot_height_all{filenum} = robot_height;
        clear td_list;
    end

% split trial data into two (based on ana_var), remove trials with missing
% markers. Store trial start and end in td_beg and td_end

    for filenum = 1:length(td_all)
        td_list = td_all{filenum};
        
        % button is > 4 V to start, then < 3.5 V in second half
        button_trial = td_list{1}.ana_var(trial_start + 5);
        trans_trial = find(button_trial < mean(button_trial),1,'first');
        trans_point = trial_start(trans_trial) - ceil(trial_start(trans_trial)-trial_end(trans_trial-1))/2;
        
        
        % split td in half based on trans_point
        fnames = fieldnames(td_list{1});
        cont_len = length(td_list{1}.ana_var);
        
        td_beg = []; td_end = [];
        for i_name = 1:numel(fnames)
            if(length(td_list{1}.(fnames{i_name})) == cont_len)
                td_beg.(fnames{i_name}) = td_list{1}.(fnames{i_name})(1:trans_point,:);
                td_end.(fnames{i_name}) = td_list{1}.(fnames{i_name})(trans_point:end,:);
            else
                td_beg.(fnames{i_name}) = td_list{1}.(fnames{i_name});
                td_end.(fnames{i_name}) = td_list{1}.(fnames{i_name});
            end
        end
        
        % put trial start and trial end into td_beg and end
        td_beg.trial_start = trial_start(1:trans_trial-1);
        td_end.trial_start = trial_start(trans_trial:end) - trans_point + 1;
        
        td_beg.trial_end = trial_end(1:trans_trial-1);
        td_end.trial_end = trial_end(trans_trial:end) - trans_point + 1;
            
        % remove trials with bad tracking
        if(isfield(td_beg,'dlc_pos'))
            for i = 1:2
                switch i
                    case 1
                        td_use = td_beg;
                    case 2
                        td_use = td_end;
                end
                % remove time points where dlc tracking is bad
                dlc_idx = [find((strcmpi(td_use.dlc_pos_names,['hand2','_x']))),...
                            find((strcmpi(td_use.dlc_pos_names,['elbow1','_x'])))];
                bad_points = any(isnan(td_use.dlc_pos(:,dlc_idx)),2) | any(isnan(td_use.dlc_vel(:,dlc_idx)),2);

                % remove trials from trial_start, trial_end where markers are
                % missing
                trial_keep = ones(numel(td_use.trial_start),1);

                for i_trial = 1:numel(td_use.trial_start)
                    if(any(bad_points(td_use.trial_start(i_trial):td_use.trial_end(i_trial)) == 1))
                        trial_keep(i_trial) = 0;
                    end
                end
                
                switch i
                    case 1
                        td_beg.trial_start = td_use.trial_start(trial_keep==1);
                        td_beg.trial_end = td_use.trial_end(trial_keep==1);
                    case 2
                        td_end.trial_start = td_use.trial_start(trial_keep==1);
                        td_end.trial_end = td_use.trial_end(trial_keep==1);
                end
            end
        end
    end
    
    td_list = {td_beg, td_end};
%% cross correlation with and without weight, only during trials

    tree_perm_use = []; max_corr = -100;
    for space_idx = 1:numel(td_list)
        % get td from td_list based on task_list
        td_use = td_list{space_idx};

        % get correlations across neurons
        exp_mask = any(1:1:numel(td_use.ana_var) >= td_use.trial_start & 1:1:numel(td_use.ana_var) <= td_use.trial_end,1);
        data_use = td_use.LeftS1_FR(exp_mask,:);
        corr_mat = corr(data_use);

        % cluster correlations  if correct task
        if(space_idx == 1)
            T = linkage(corr_mat);
            f_to_close = figure();
            [~,~,tree_perm_use] = dendrogram(T);
            close(f_to_close);
        end
        % get max corr for color limits later
        for i = 1:size(corr_mat,1) corr_mat(i,i) = nan; end % set diagonal as nan
        max_corr = max(max_corr,max(max(corr_mat)));
    end

    corr_mat_use = {};
    for space_idx = 1:numel(td_list)
        % get td from td_list based on task_list
        td_use = td_list{space_idx};
        % plot correlations
        f=figure();
        exp_mask = any(1:1:numel(td_use.ana_var) > td_use.trial_start & 1:1:numel(td_use.ana_var) < td_use.trial_end,1);
        data_use = td_use.LeftS1_FR(exp_mask,:);
        corr_mat = corr(data_use);
        corr_mat_use{space_idx} = corr(data_use(:,tree_perm_use));
        for i = 1:size(corr_mat_use{space_idx},1) corr_mat_use{space_idx}(i,i) = nan; end % set diagonal as nan
        imagesc(corr_mat_use{space_idx});
        b=colorbar;
        caxis(max_corr*[-1,1]);
        cmap = colorcet('D1'); colormap(cmap);
    end
        

    f=figure();
    corr_mat_diff = corr_mat_use{2} - corr_mat_use{1};
    for i = 1:size(corr_mat_use{space_idx},1) corr_mat_use{space_idx}(i,i) = nan; end % set diagonal as nan
    imagesc(corr_mat_diff);
    b=colorbar;
    caxis(max(max(corr_mat_diff,[],'omitnan'),[],'omitnan')*[-1,1]);
    cmap = colorcet('D1'); colormap(cmap);

        
%% Loop through files to cross-validate encoders
    
    % build one model for the entire data. Include term for whether the
    % weight was on or not
    
    td_use = [];
    
    for i_td = 1:numel(td_list)
        % get keep_mask
        keep_mask = any(1:1:numel(td_list{i_td}.ana_var) >= td_list{i_td}.trial_start & 1:1:numel(td_list{i_td}.ana_var) <= td_list{i_td}.trial_end,1);
        % put relevant data into td_all
        if(i_td==1)
            td_beg_size = sum(keep_mask);

            % data
            td_use.LeftS1_FR = td_list{i_td}.LeftS1_FR(keep_mask==1,:);
            td_use.dlc_pos = td_list{i_td}.dlc_pos(keep_mask==1,:);
            td_use.dlc_vel = td_list{i_td}.dlc_vel(keep_mask==1,:);
            td_use.weight_var = ones(sum(keep_mask),1);
            
            % meta info
            td_use.dlc_pos_names = td_list{i_td}.dlc_pos_names;
            td_use.LeftS1_unit_guide = td_list{i_td}.LeftS1_unit_guide;
            td_use.bin_size = td_list{i_td}.bin_size'
        else
            % append data
            td_use.LeftS1_FR = [td_use.LeftS1_FR;td_list{i_td}.LeftS1_FR(keep_mask==1,:)];
            td_use.dlc_pos = [td_use.dlc_pos;td_list{i_td}.dlc_pos(keep_mask==1,:)];
            td_use.dlc_vel = [td_use.dlc_vel;td_list{i_td}.dlc_vel(keep_mask==1,:)];
            td_use.weight_var = [td_use.weight_var; zeros(sum(keep_mask),1)];
        end
    end
    
    % build glm with hand+elbow pos+vel as predictors. Include weight term
    % and bias as predictors as well
    % indices for cartesian hand coordinates
    %
    % smooth firing rates
    td_use = smoothSignals(td_use,struct('signals','LeftS1_FR'));
    
    markername = 'hand2';
    [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_use.dlc_pos_names);
    assert(all(point_exists),'Hand marker does not exist?')

    markername = 'elbow1';
    [point_exists,marker_elbow_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_use.dlc_pos_names);
    assert(all(point_exists),'Elbow marker does not exist?')

    x = [td_use.dlc_pos(:,[marker_hand_idx marker_elbow_idx]), ...
        td_use.dlc_vel(:,[marker_hand_idx marker_elbow_idx]),...
            td_use.weight_var(:)];
        
    y = td_use.LeftS1_FR(:,:);
    
    glm_distribution = 'poisson';
    b = zeros(size(x,2)+1,size(y,2));
    yfit = zeros(size(y));
    
    for iVar = 1:size(y,2)
        [b(:,iVar),~,s_temp] = glmfit(x,y(:,iVar),glm_distribution);
        yfit(:,iVar) = exp([ones(size(x,1),1), x]*b(:,iVar)); 
    end
    
    
    