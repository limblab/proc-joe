%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% freereaching_analysis - 
%%      script to perform the free reach vs constrained analyses.
%%      Note: this script is essentially split into two parts:
%%      The first part fits the encoder models using cross-validation
%%      while the second part plots the results. This segregation is
%%      for ease of use, as the cross-validation can take a significant
%%      amount of time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up meta info and load trial data
    if ispc
        folderpath = 'D:\Lab\Data\DLC_videos\Han_20201204_rwFreeReach\neural-data';
    else
        folderpath = '/data/raeed/project-data/limblab/s1-kinematics';
    end
    
    % load data
    file_info = dir(fullfile(folderpath,'*constrainedExp_td*'));
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
    monkey_names = {'H'};
    included_models = {'ext','handelbow'}; % models to calculate encoders for
    models_to_plot = {'ext','handelbow'}; % main models of the paper
    not_plot_models = setdiff(included_models,models_to_plot);
    bin_size = 0.05; % s
    
    % colors for pm, dl conditions and sessions
    cond_colors = [...
        231,138,195;...
        166,216,84]/255;
    session_colors = [...
        102,194,165;...
        252,141,98;...
        141,160,203]/255;

%% Loop through trial data files to clean up
    td_all = {};
    for filenum = 1:length(filenames)
        % Load data
        load(fullfile(file_info(filenum).folder,file_info(filenum).name));
        
        for i_td = 1:numel(td_list) % each entry in td_list is a different trial_data for a different experiment (free reach vs 2D random walk for example)            
            % resample trial data to appropriate bin size
            if(td_list{i_td}.bin_size <= bin_size)
                td_list{i_td} = binTD(td_list{i_td},bin_size/td_list{i_td}.bin_size);
            else
                warning('td bin size is larger than desired bin size');
            end
            % get marker velocity
            td_list{i_td} = getDifferential(td_list{i_td},struct('signals','dlc_pos','alias','dlc_vel'));
        
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
            
            % remove time points where dlc tracking is bad
            bad_points = any(isnan(td_list{i_td}.dlc_pos),2) | any(isnan(td_list{i_td}.dlc_vel),2);
            td_names = fieldnames(td_list{i_td});
            for i_name = 1:numel(td_names)
                if(size(bad_points,1) == size(td_list{i_td}.(td_names{i_name}),1))
                    td_list{i_td}.(td_names{i_name})(bad_points,:) = [];
                end
            end
            fprintf('Removed %d percent of trials because of missing markers\n',sum(bad_points)/numel(bad_points))
            
            % remove trials where monkey's arm position is out of the
            % "workspace"
             
        end
        td_all{filenum} = td_list;
    end

    
    
%% compute cross correlation between firing rate and first pc of velocity data
    markername = 'hand3';
    for filenum = 1:numel(td_all)
        td_list = td_all{filenum};
        best_neural_lag = nan(numel(td_list),numel(td_list{1}.([arrayname,'_ts'])));
        best_neural_corr = nan(size(best_neural_lag));
        
        best_handle_lag = nan(numel(td_list),1);
        best_handle_corr = nan(numel(td_list),1);
        
        for i_td = 1%:numel(td_list)
            % get first pc of vel data for a single marker (something on
            % the hand)
            dlc_idx = find((strcmpi(td_list{1}.dlc_names,[markername,'_x'])));
            vel_data = td_list{i_td}.dlc_vel(:,[dlc_idx,dlc_idx+1,dlc_idx+2]);
            [coeff,score,latent] = pca(vel_data);
            dlc_pc1_data = score(:,1);
            
            % get correlation bewteen neural data and first pc
            for i_unit = 1:size(td_list{i_td}.([arrayname,'_unit_guide']),1)
                [r,lags] = xcorr(td_list{i_td}.([arrayname,'_spikes'])(:,i_unit),dlc_pc1_data,200,'coeff');
                [~,best_idx] = max(abs(r));
                best_neural_lag(i_td,i_unit) = lags(best_idx);
                best_neural_corr(i_td,i_unit) = r(best_idx);
            end
            
            % get correlation between marker and handle position (first pc
            % for both)
            if(strcmpi(task_list{i_td},'RT'))
                vel_data = td_list{i_td}.vel;
                [coeff,score,latent] = pca(vel_data);
                handle_pc1_data = score(:,1);
                [r_handle, lags_handle] = xcorr(dlc_pc1_data,handle_pc1_data,200,'coeff');
                [~,best_idx] = max(abs(r_handle));
                best_handle_lag = lags_handle(best_idx);
                best_handle_corr = r_handle(best_idx);
            end
            
        end
        
    end

    
%% Loop through files to cross-validate encoders
%     fprintf('Starting analysis of %d files. Warning: cross-validation will take a long time\n',length(trial_data_cell))
%     encoderResults_cell = cell(length(trial_data_cell),1);
%     for filenum = 1:length(trial_data_cell)
%         % Load data
%         td = trial_data_cell{filenum};
%     
%         % bin data at 50ms
%         td = binTD(td,0.05/td(1).bin_size);
%     
%         %% Get encoding models
%         encoderResults_cell{filenum} = mwEncoders(td,struct(...
%             'model_aliases',{included_models},...
%             'arrayname',arrayname,...
%             'num_tuning_bins',16,...
%             'crossval_lookup',[],...
%             'get_tuning_curves',true,...
%             'num_repeats',20,...
%             'num_folds',5));
%     end
% 
% %% save encoder results
%     if savefile
%         fprintf('Saving encoder results files...\n')
%         save(fullfile(savedir,savename),'encoderResults_cell')
%     end
