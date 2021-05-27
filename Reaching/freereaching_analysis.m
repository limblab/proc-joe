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
    clear; close all; clc;
    if ispc
        %base_folderpath = 'C:\Users\dongq\DeepLabCut\Crackle-Qiwei-2020-12-03\';
        base_folderpath = 'C:\Users\dongq\DeepLabCut\Crackle-Qiwei-2020-12-15\';
        %base_folderpath = 'C:\Users\dongq\DeepLabCut\Crackle-Qiwei-2020-12-16\';
        %base_folderpath = 'C:\Users\dongq\DeepLabCut\Han_20201203_rwFreeReach\'; %PROBLEM: RT3D task does not have trial_idx
        %base_folderpath = 'C:\Users\dongq\DeepLabCut\Han_20201204_rwFreeReach\'; %PROBLEM: RT3D task does not have trial_idx
        %base_folderpath = 'C:\Users\dongq\DeepLabCut\Han_20201217_rwFreeReach\';
        folderpath = strcat(base_folderpath, 'neural-data\');
        %folderpath = 'C:\Users\dongq\DeepLabCut\!DLC folders not to be used for thesis but can be used later\Han_20201204_rwFreeReach\neural-data\';
        %folderpath = 'C:\Users\dongq\DeepLabCut\Han_20201204_rwFreeReach\neural-data\';
        %folderpath = 'D:\DLCdata\Han_20201204_rwFreeReach\neural-data\';
        %folderpath = 'C:\Users\dongq\DeepLabCut\Han_20201203_rwFreeReach\neural-data\';
        %folderpath = 'C:\Users\dongq\DeepLabCut\Han_20201217_rwFreeReach\neural-data\';
        %folderpath = 'C:\Users\dongq\DeepLabCut\!DLC folders not to be used for thesis but can be used later\Crackle-Qiwei-2020-12-03\neural-data\';

        %folderpath = 'C:\Users\dongq\DeepLabCut\Crackle-Qiwei-2020-12-15\neural-data\';
        %folderpath = 'C:\Users\dongq\DeepLabCut\Crackle-Qiwei-2020-12-16\neural-data\';
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
    %monkey_names = {'Han'};
    monkey_names = {'Crackle'};
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
                td_list{i_td} = setOriginAsShoulder(td_list{i_td},1);
                % get marker velocity
                td_list{i_td} = getDifferential(td_list{i_td},struct('signals','dlc_pos','alias','dlc_vel'));
                % remove time points where dlc tracking is bad
                dlc_idx = [find((strcmpi(td_list{i_td}.dlc_pos_names,['hand2','x']))),...
                            find((strcmpi(td_list{i_td}.dlc_pos_names,['hand2','_y']))),...
                            find((strcmpi(td_list{i_td}.dlc_pos_names,['elbow1','_x']))),...
                            find((strcmpi(td_list{i_td}.dlc_pos_names,['elbow1','_y'])))];
                bad_points = any(isnan(td_list{i_td}.dlc_pos(:,dlc_idx)),2) | any(isnan(td_list{i_td}.dlc_vel(:,dlc_idx)),2);
                td_names = fieldnames(td_list{i_td});
                for i_name = 1:numel(td_names)
                    if(size(bad_points,1) == size(td_list{i_td}.(td_names{i_name}),1))
                        td_list{i_td}.(td_names{i_name})(bad_points,:) = [];
                    end
                end
                fprintf('Removed %.2f%% percent of trials because of missing markers\n',sum(bad_points)/numel(bad_points)*100)
            end
            % remove trials where monkey's arm position is out of the
            % "workspace"
            %td_list{i_td}
            %task_list{i_td}
            td_list{i_td} = getExperimentPhase(td_list{i_td},task_list{i_td});
% 
% 
%             % get robot height (z data of a hand marker). This is
%             % meaningless during free reach...
%             markername = 'hand3';
%             dlc_idx = find((strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_z'])));
%             robot_height(end+1) = mean(td_list{i_td}.dlc_pos(:,dlc_idx),'omitnan');
        end
        td_all{filenum} = td_list;
        task_list_all{filenum} = task_list;
        robot_height_all{filenum} = robot_height;
        clear td_list;
    end

    
    
%% Loop through files to make kinematic figures
    fprintf('Starting kinematic analysis of %d files.',length(td_all))
    kin_input_data = [];
    for filenum = 1:length(td_all)
        td_list = td_all{filenum};
        task_list = task_list_all{filenum};
        
        % get kinematic data and make plots
        kin_data = reachingKinematics(td_list,task_list,kin_input_data);
    end
    
%% get pca of muscle lengths for predictions
    for filenum = 1:length(td_all)
        td_list = td_all{filenum};
        for i_td = 1:numel(td_list)
            PCAparams = struct('signals',{{'opensim',find(contains(td_list{i_td}.opensim_names,'_len') & ~contains(td_list{i_td}.opensim_names,'tricep_lat'))}},...
                                    'do_plot',true);
    %                 PCAparams = struct('signals','muscle_len', 'do_plot',false);
            [td_temp,~] = dimReduce(td_list{i_td},PCAparams);
            td_list{i_td}.musc_len_pca = td_temp.opensim_pca;
            % get velocity PCA
            % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
            PCAparams_vel = struct('signals',{{'opensim',find(contains(td_list{i_td}.opensim_names,'_muscVel') & ~contains(td_list{i_td}.opensim_names,'tricep_lat'))}},...
                                'do_plot',true);
            [td_temp,~] = dimReduce(td_list{i_td},PCAparams_vel);
            td_list{i_td}.musc_vel_pca = td_temp.opensim_pca;
            clear td_temp;
        end
        td_all{filenum} = td_list;
    end

    
%% Loop through files to cross-validate encoders
    fprintf('Starting analysis of %d files. Warning: cross-validation will take a long time\n',length(td_all))
    encoderResults_cell = cell(length(td_all),1);
    for filenum = 1:length(td_all)
        % Load data
        td_list = td_all{filenum};
        task_list = task_list_all{filenum};
        % bin data at 50ms
        for i_td = 1:numel(td_list)
            td_list{i_td} = binTD(td_list{i_td},0.05/td_list{i_td}(1).bin_size);
        end
        
        % Get encoding models
        encoderResults_cell{filenum} = reachingEncoders(td_list,task_list,robot_height_all{filenum},struct(...
            'model_aliases',{included_models},...
            'arrayname',arrayname,...
            'num_tuning_bins',16,...
            'crossval_lookup',[],...
            'get_tuning_curves',true,...
            'num_repeats',5,... % Raeed used 20 repeats.
            'num_folds',5));
    end

%% save encoder results
    if savefile
        fprintf('Saving encoder results files...\n')
        save(fullfile(savedir,savename),'encoderResults_cell')
    end
    
%% load encoder results (assuming a break)
    if savefile
        fprintf('Loading encoder results files...\n')
        load(fullfile(savedir,savename))
    end

%% Compile important information over all files
    [model_eval,model_tuning,tuning_corr,shift_vaf,tuned_neurons] = deal(cell(length(monkey_names),size(session_colors,1)));
    session_ctr = zeros(length(monkey_names),1);
    fileclock = tic;
    fprintf('Started loading files...\n')
    for filenum = 1:length(encoderResults_cell)
        % load data
        encoderResults = encoderResults_cell{filenum};

        % classify monkey and session number
        monkey_idx = find(strcmpi(encoderResults.crossEval.monkey{1},monkey_names));
        session_ctr(monkey_idx) = session_ctr(monkey_idx) + 1;

        % We already have evaluation table in crossEval... just extract the models we want
        model_eval{monkey_idx,session_ctr(monkey_idx)} = encoderResults.crossEval(:,contains(encoderResults.crossEval.Properties.VariableDescriptions,'meta'));
        model_eval_cell = cell(1,length(included_models));
        [space_eval_cell,space_eval_within_cell,space_eval_across_cell] = deal(cell(2,length(included_models)));
        for modelnum = 1:length(included_models)
            model_eval_cell{modelnum} = table(encoderResults.crossEval.(sprintf('glm_%s_model_eval',included_models{modelnum})),...
                'VariableNames',strcat(included_models(modelnum),'_eval'));
            model_eval_cell{modelnum}.Properties.VariableDescriptions = {'linear'};
            for spacenum = 1:2
                space_eval_cell{spacenum,modelnum} = table(encoderResults.crossEval.(sprintf('glm_%s_model_space%d_eval',included_models{modelnum},spacenum)),...
                    'VariableNames',{sprintf('%s_space%d_eval',included_models{modelnum},spacenum)});
                space_eval_cell{spacenum,modelnum}.Properties.VariableDescriptions = {'linear'};
                % because some old files don't have this...
                try
                    space_eval_within_cell{spacenum,modelnum} = table(encoderResults.crossEval.(sprintf('glm_%s_model_space%d_within_eval',included_models{modelnum},spacenum)),...
                        'VariableNames',{sprintf('%s_space%d_within_eval',included_models{modelnum},spacenum)});
                    space_eval_within_cell{spacenum,modelnum}.Properties.VariableDescriptions = {'linear'};
                catch ME
                    warning('Within space predictions are not available. Eval table is not completely filled out')
                end
                
                try
                    space_eval_across_cell{spacenum,modelnum} = table(encoderResults.crossEval.(sprintf('glm_%s_model_space%d_across_eval',included_models{modelnum},spacenum)),...
                        'VariableNames',{sprintf('%s_space%d_across_eval',included_models{modelnum},spacenum)});
                    space_eval_across_cell{spacenum,modelnum}.Properties.VariableDescriptions = {'linear'};
                catch ME
                    warning('Across space predictions are not available. Eval table is not completely filled out')
                end
                
            end
        end
        model_eval{monkey_idx,session_ctr(monkey_idx)} = horzcat(...
            model_eval{monkey_idx,session_ctr(monkey_idx)},...
            model_eval_cell{:},...
            space_eval_cell{:},...
            space_eval_within_cell{:},...
            space_eval_across_cell{:});

        % We already have tuning table in crossTuning... just extract the models we want
        model_tuning{monkey_idx,session_ctr(monkey_idx)} = encoderResults.crossTuning(:,...
            contains(encoderResults.crossTuning.Properties.VariableDescriptions,'meta') |...
            strcmpi(encoderResults.crossTuning.Properties.VariableNames,'bins'));
        model_tuning_cell = cell(1,length(included_models)+1);
        for modelnum = 1:length(included_models)
            model_tuning_cell{modelnum} = table(...
                encoderResults.crossTuning.(sprintf('glm_%s_model_dlc_vel_handxyCurve',included_models{modelnum})),...
                encoderResults.crossTuning.(sprintf('glm_%s_model_dlc_vel_handxyPD',included_models{modelnum})),...
                'VariableNames',strcat(included_models(modelnum),{'_velCurve','_velPD'}));
            model_tuning_cell{modelnum}.Properties.VariableDescriptions = {'linear','circular'};
        end
        model_tuning_cell{end} = table(...
            encoderResults.crossTuning.([arrayname,'_FR_dlc_vel_handxyCurve']),...
            encoderResults.crossTuning.([arrayname,'_FR_dlc_vel_handxyPD']),...
            'VariableNames',strcat('S1_FR',{'_velCurve','_velPD'}));
        model_tuning_cell{end}.Properties.VariableDescriptions = {'linear','circular'};
        % put it together
        model_tuning{monkey_idx,session_ctr(monkey_idx)} = horzcat(...
            model_tuning{monkey_idx,session_ctr(monkey_idx)},...
            model_tuning_cell{:});

        % Get tuning curve correlation table
        tuning_corr{monkey_idx,session_ctr(monkey_idx)} = calculateEncoderTuningCorr(...
            encoderResults,struct('model_aliases',{included_models},'neural_signal',[arrayname,'_FR']));

        % get tuned neurons
        if isfield(encoderResults,'tunedNeurons')
            % Get PD shift error table
            shift_vaf{monkey_idx,session_ctr(monkey_idx)} = calculateEncoderPDShiftVAF(...
                encoderResults,struct('model_aliases',{included_models}));
            tuned_neurons{monkey_idx,session_ctr(monkey_idx)} = encoderResults.tunedNeurons;
        else
            warning('No tuned neurons field found!')
        end
        
        
        % output a counter
        fprintf('Processed file %d of %d at time %f\n',filenum,length(encoderResults_cell),toc(fileclock))
    end

%% plot random "trials" and predictions from one of the models - real FR, model predict FR

    unit_idx = 10;
    td_idx = 2; % pick which trial data (task) to use

    window_plot = [-2,2]; % s
    num_trials_plot = 3;
    bin_size = td_list{td_idx}.bin_size;
    % pick random time point, make sure data in window are consecutive
    % based on trial_idx, if yes go on and plot
    
    td_list_smooth = smoothSignals(td_list{td_idx},'LeftS1_FR');

    time_idx = []; window_idx = [];
    for i_trial = 1:num_trials_plot
        is_consec = 0;
        test_ctr = 0;
        while ~is_consec && test_ctr < 100
            time_idx(i_trial) = (rand()*(numel(td_list{td_idx}.trial_idx) - floor(diff(window_plot)*2/bin_size))) + floor(diff(window_plot)/bin_size);
            window_idx(i_trial,:) = time_idx(i_trial)+floor(window_plot./bin_size);
            trial_idx = td_list{td_idx}.trial_idx(window_idx(i_trial,1):window_idx(i_trial,2));

            if(max(diff(trial_idx)) == 1)
                is_consec = 1;
            end

            test_ctr = test_ctr + 1;
        end
    end
    
    %figure('Position',[145 558 1095 420]); hold on;
    figure('Position',[145 558 1840 300]); hold on;
    % plot FR for unit
    x_data = (0:1:diff(window_idx(1,:)))*bin_size;
    encoderResults = encoderResults_cell{1};
    repeatnum = 1;
    % predict FR for unit in time window, then plot for each model
    for mdlnum = 1:numel(models_to_plot)
        td_list_smooth = getModel(td_list_smooth,encoderResults.glm_info{repeatnum,mdlnum});
    end
    for i_trial = 1:num_trials_plot
        %i_trial
        ax_list(i_trial) = subplot(1,num_trials_plot,i_trial); hold on;
        plot(x_data,td_list_smooth.LeftS1_FR(window_idx(i_trial,1):window_idx(i_trial,2),unit_idx),'k','linewidth',2)

        for mdlnum = 1:numel(models_to_plot)
            glm_fieldname = ['glm_',models_to_plot{mdlnum},'_model'];
            plot(x_data,td_list_smooth.(glm_fieldname)(window_idx(i_trial,1):window_idx(i_trial,2),unit_idx),'color',getColorFromList(1,mdlnum-1),'linewidth',2)
        end
        
        formatForLee(gcf);
        set(gca,'fontsize',14);
        
    end
    
    linkaxes(ax_list,'xy');
    subplot(1,num_trials_plot,1); 
    l=legend('Actual','Hand-only','Whole-arm');
    set(l,'box','off');
    
    xlabel('Time (s)');
    ylabel('FR (Hz)');
    title(td_list{td_idx}.monkey + " date " + td_list{td_idx}.date + " task " + td_list{td_idx}.task + " neuron " + string(unit_idx));
   

    
%% %% plot random "trials" and predictions from one of the models - real FR, model predict FR ON MULTIPLE NEURONS

    unit_idx = 10;
    %unit_idx_list = [10,11,12]; %RT3D
    %For Han-1204
    %unit_idx_list = [13,19,25]; %RT3D 13,16,19,22,25,31
    %For Crackle-1215
    unit_idx_list = [4,14,18];
    %For Crackle-1216
    %unit_idx_list = [10,11,14];
    
    %unit_idx_list = [1,2,3]; % 1,2 not responding
    %unit_idx_list = [4,5,6]; %6 not responding
    %unit_idx_list = [7,8,9]; %8,9 not responding (I actually want to say
    %"modulating")
    %unit_idx_list = [10,11,12];
    %unit_idx_list = [13,14,15]; %14, 15 not responding
    %unit_idx_list = [16,17,18]; %17, 18 kinda not responding
    %unit_idx_list = [19,20,21];
    %unit_idx_list = [22,23,24]; %23,24 not responding
    %unit_idx_list = [25,26,27]; %27 not responding
    %unit_idx_list = [28,29,30]; %30 not responding
    %unit_idx_list = [31];
    %unit_idx_list = [12,25,31]; %RT2D
    num_rows = length(unit_idx_list) + 1 %last row for speed 
    td_idx = 1; % pick which trial data (task) to use, 1 for RT2D 2 for RT3D

    window_plot = [-2,2]; % s
    num_trials_plot = 3;
    num_neurons = length(unit_idx_list)
    bin_size = td_list{td_idx}.bin_size;
    % pick random time point, make sure data in window are consecutive
    % based on trial_idx, if yes go on and plot
    
    td_list_smooth = smoothSignals(td_list{td_idx},'LeftS1_FR');

    %randomly choose a time point
%     time_idx = []; window_idx = [];
%     for i_trial = 1:num_trials_plot
%         is_consec = 0;
%         test_ctr = 0;
%         while ~is_consec && test_ctr < 100
%             time_idx(i_trial) = (rand()*(numel(td_list{td_idx}.trial_idx) - floor(diff(window_plot)*2/bin_size))) + floor(diff(window_plot)/bin_size);
%             window_idx(i_trial,:) = time_idx(i_trial)+floor(window_plot./bin_size);
%             trial_idx = td_list{td_idx}.trial_idx(window_idx(i_trial,1):window_idx(i_trial,2));
% 
%             if(max(diff(trial_idx)) == 1)
%                 is_consec = 1;
%             end
%             
%             
% 
%             test_ctr = test_ctr + 1;
%         end
%     end

    window_length = ( window_plot(2) - window_plot(1) )/bin_size
%     RT3D_section1 = 8600;
%     RT3D_section2 = 1230;
%     RT3D_section3 = 14000;
%     window_idx = [RT3D_section1,RT3D_section1+window_length; 
%         RT3D_section2,RT3D_section2+window_length; 
%         RT3D_section3,RT3D_section3+window_length];
    
    RT2D_section1 = 8000;
    RT2D_section2 = 6980;
    RT2D_section3 = 2400;
    
    %Try for Crackle-1216
%     RT2D_section1 = 7000;
%     RT2D_section2 = 500;
%     RT2D_section3 = 4020;
    
    window_idx = [RT2D_section1,RT2D_section1+window_length; 
        RT2D_section2,RT2D_section2+window_length; 
        RT2D_section3,RT2D_section3+window_length];

    %figure('Position',[145 558 1095 420]); hold on;
    figure('Position',[10 50 1840 700]); hold on;
    % plot FR for unit
    x_data = (0:1:diff(window_idx(1,:)))*bin_size;
    encoderResults = encoderResults_cell{1};
    repeatnum = 1;
    % predict FR for unit in time window, then plot for each model
    for mdlnum = 1:numel(models_to_plot)
        td_list_smooth = getModel(td_list_smooth,encoderResults.glm_info{repeatnum,mdlnum});
    end
    for i_neuron = 1:num_rows
        for j_trial = 1:num_trials_plot
            %i_trial
            current_plot_num = (i_neuron - 1) * num_trials_plot + j_trial
            ax_list(j_trial) = subplot(num_rows,num_trials_plot,current_plot_num); hold on;
            
            %plot the actual firing rate data
            if i_neuron ~= num_rows
                actual_y_data = td_list_smooth.LeftS1_FR(window_idx(j_trial,1):window_idx(j_trial,2),unit_idx_list(i_neuron));
                plot(x_data,actual_y_data,'k','linewidth',2)
            end
            
            if j_trial == 1 & i_neuron ~= num_rows%the first column
                ylabel("neuron " + string(unit_idx_list(i_neuron)));
            end
            
            %Plot speed
            if i_neuron == num_rows %if last row
                dlc_section_pos = td_list_smooth.dlc_pos(window_idx(j_trial,1):window_idx(j_trial,2),19:21);
                dlc_section_velocity = diff(dlc_section_pos);
                dlc_section_speed = sqrt(dlc_section_velocity(:,1).^2 + dlc_section_velocity(:,2).^2 + dlc_section_velocity(:,3).^2)./td_list_smooth.bin_size/100;
                dlc_section_speed(length(dlc_section_speed)+1,:) = dlc_section_speed(length(dlc_section_speed),:);
                dlc_section_acc = diff(dlc_section_speed);
                dlc_section_acc(length(dlc_section_acc)+1,:) = dlc_section_acc(length(dlc_section_acc),:);
                %plot(x_data,dlc_section_speed*10,'color',getColorFromList(1,mdlnum),'linewidth',2)
                plot(x_data,dlc_section_speed,'color',getColorFromList(1,mdlnum),'linewidth',2);
                %plot(x_data,dlc_section_acc,'linewidth',2);
                %ylim([0,2.5])
                %plot(x_data,td_list_smooth.vel(window_idx(j_trial,1):window_idx(j_trial,2),1));
                %plot(x_data,td_list_smooth.vel(window_idx(j_trial,1):window_idx(j_trial,2),2));
                %plot(x_data,td_list_smooth.acc(window_idx(j_trial,1):window_idx(j_trial,2),1));
                %plot(x_data,td_list_smooth.acc(window_idx(j_trial,1):window_idx(j_trial,2),2));
                if j_trial == 1
                    ylabel("speed (in m/s)")
                end
            end
            
            
            %try to calculate pR2 here in this section where I plot the
            %neural data
            %plot estimated neural data
            for mdlnum = 1:numel(models_to_plot)
                 if i_neuron ~= num_rows %if not the last row
                %if i_neuron == num_rows % if is the last row
                    glm_fieldname = ['glm_',models_to_plot{mdlnum},'_model'];
                    y_data = td_list_smooth.(glm_fieldname)(window_idx(j_trial,1):window_idx(j_trial,2),unit_idx_list(i_neuron));
                    plot(x_data,y_data,'color',getColorFromList(1,mdlnum-1),'linewidth',2)
                    %ylim([0,100])
%                     %calculate pR2
%                     [~, dev] = glmfit(x_data, y_data, 'binomial');
%                     [~, devModel] = glmfit(x_data, y_data, 'binomial'); 
%                     [~, devNull]  = glmfit(x_data, y_data(randperm(length(y_data))), 'binomial');
%                     R2 = (DevNull-DevModel)/(DevNull)
                end
%                 if i_neuron == num_rows % if is the last row
%                     temp_1d_actual_y_data = reshape(actual_y_data.',1,[]);
%                     temp_1d_y_data = reshape(y_data.',1,[]);
%                     temp_pR2 = pseudoR2(temp_1d_actual_y_data, temp_1d_y_data, mean(actual_y_data));
%                     
%                     %temp_1d_actual_y_data(end+1) = temp_1d_actual_y_data(end);
%                     %temp_1d_y_data(end+1) = temp_1d_y_data(end);
%                     
%                     temp_pR2_single_frame = [];
%                     temp_1d_actual_y_data(1)
%                     temp_1d_y_data(1)
%                     for i = 1:numel(temp_1d_actual_y_data)
%                         temp_pR2_single_frame(i) = pseudoR2(temp_1d_actual_y_data(i), temp_1d_y_data(i), temp_1d_actual_y_data(i));
%                         %temp_pR2_single_frame(i) = pseudoR2(temp_1d_actual_y_data(i), temp_1d_y_data(i), mean(temp_1d_actual_y_data));
%                     end
%                     %temp_pR2_single_frame(1)
%                     plot(x_data, temp_pR2_single_frame,'color',getColorFromList(1,mdlnum-1),'linewidth',2)
%                     
%                     
%                 end
            end

            formatForLee(gcf);
            set(gca,'fontsize',14);
            
            %don't plot the x and y ticks if not the first column
            if j_trial ~= 1
                set(gca, 'YTick', [])
            end
            if i_neuron ~= num_rows
                set(gca,'XTick',[])
            end
            
        end
    end
    %linkaxes(ax_list,'xy');
    subplot(num_rows,num_trials_plot,1); 
    %l=legend('Actual','Hand-only','Whole-arm');
    %set(l,'box','off');
    
    xlabel('Time (s)');
    ylabel('FR (Hz)');
    %title(td_list{td_idx}.monkey + " date " + td_list{td_idx}.date + " task " + td_list{td_idx}.task + " neuron " + string(unit_idx));
    title(td_list{td_idx}.monkey + " date " + td_list{td_idx}.date + " task " + td_list{td_idx}.task);
        
    
    subplot(num_rows,num_trials_plot,2); 
    ylabel("neuron " + unit_idx_list(1))
    
    
    
%% Plot spike raster plots with hand speed
% task = 2 %1 for RT2D, 2 for RT3D
% %hand 2 is usually on column 19~21
% [num_of_bins,num_of_neurons] = size(td_list{task}.LeftS1_FR);
% 
% %neurons_to_use = [10,11,12,19,20,21,22,26,28,29,31];
% neurons_to_use = 1:num_of_neurons;
% 
% num_of_neurons = length(neurons_to_use);
% 
% example_section_start = 5924;
% example_section_length = 50 %num of 0.05s bins
% example_section_start = randi([0 num_of_bins-example_section_length],1);
% example_section_end = example_section_start + example_section_length;
% 
% x_data_section = zeros((example_section_end-example_section_start)+1,num_of_neurons);
% for i = 1:num_of_neurons
%     x_data_section(:,i) = (0:1:(example_section_end-example_section_start))*bin_size;%transpose((0:1:(example_section_end-example_section_start))*bin_size)
% end
% %example_raster_section = td_list{task}.LeftS1_spikes(example_section_start:example_section_end,:);
% example_raster_section = td_list{task}.LeftS1_spikes(example_section_start:example_section_end,neurons_to_use);
% example_raster_section(example_raster_section>1) = 1;
% example_raster_timepoints = x_data_section.*example_raster_section;
% example_raster_section_hist = td_list{task}.LeftS1_spikes(example_section_start:example_section_end,neurons_to_use);
% example_raster_section_hist = sum(example_raster_section_hist);
% 
% 
% % Change the datastructure to what plot_raster.m wants
% % x array with the data of all the time points that a spike occured
% % y array with the data of which neuron this spike belongs to 
% example_spike_times = [];
% example_spike_times_neuron = [];
% for i = 1:num_of_neurons
%     for j = 1:length(example_raster_timepoints)
%         if example_raster_timepoints(j,i) ~= 0
%             example_spike_times(end+1) = example_raster_timepoints(j,i);
%             example_spike_times_neuron(end+1) = i;
%         end
%     end
% end
% %x_data_section = transpose((0:1:(example_section_end-example_section_start),0:  32)*bin_size);
% figure
% example_kinematics_section_pos = td_list{task}.dlc_pos(example_section_start:example_section_end,19:21)
% example_kinematics_section_vel = diff(example_kinematics_section_pos) %difference in position for each 0.05s bin
% example_kinematics_section_spd = sqrt(example_kinematics_section_vel(:,1).^2 + example_kinematics_section_vel(:,2).^2 + example_kinematics_section_vel(:,3).^2)./td_list_smooth.bin_size/100; %now in m/s
% example_kinematics_section_spd(length(example_kinematics_section_spd)+1,:) = example_kinematics_section_spd(length(example_kinematics_section_spd),:);
% example_kinematics_section_spd(example_kinematics_section_spd>2) = 2
% hold()
% plot_raster(example_spike_times,example_spike_times_neuron,1,[0.8500 0.3250 0.0980],1)
% 
% %barh(example_raster_section_hist);
% %speed = plot(x_data_section,example_kinematics_section_spd*5,'LineWidth',2)
% %speed = plot(x_data_section,example_kinematics_section_spd*5,'Color',[0.4940 0.1840 0.5560],'LineWidth',2)
% ylabel("each neuron")
% %xlim([0,80])
% xlabel("number of spikes in this section")
% %alpha(speed,0.5)
                    
    
    
%% Get pR2 pairwise comparisons for model pairs and all neurons
    % find winners of pR2
    pr2_winners = cell(length(monkey_names),size(session_colors,1));
    for monkeynum = 1:length(monkey_names)
        for sessionnum = 1:session_ctr(monkeynum)
            [pr2_winners{monkeynum,sessionnum},model_pairs] = compareEncoderMetrics(...
                    model_eval{monkeynum,sessionnum},struct(...
                        'bonferroni_correction',6,...
                        'models',{models_to_plot},...
                        'model_pairs',{{'ext','handelbow'}},...
                        'postfix','_eval'));
        end
    end

    % figure out how many neurons the hand models could beat either of the whole-arm models
    all_pr2_winners = horzcat(pr2_winners{:});
    ext_winners = sum(strcmpi(all_pr2_winners,'ext'),2);
    handelbow_winners = sum(strcmpi(all_pr2_winners,'handelbow'),2);
    fprintf('pR2 winners -- hand-only: %d, whole-arm: %d\n',ext_winners,handelbow_winners)

    % compare two models trained and tested within each task. Subplot for each task, each model is
    % an axis
    sessionnum = 1;
    spacenames = {'RT3D','RT2D'};
    legend_data = [];

    figure
    for monkeynum = 1:length(monkey_names)
        for spacenum = 1:2 % 1 = 3D, 2 = 2D
            % set subplot
            subplot(1,length(monkey_names),monkeynum)
            plot([-1 1],[-1 1],'k--','linewidth',0.5)
            hold on
            plot([0 0],[-1 1],'k-','linewidth',0.5)
            plot([-1 1],[0 0],'k-','linewidth',0.5)
            
            avg_pR2 = neuronAverage(model_eval{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
            % scatter filled circles if there's a winner, empty circles if not
            no_winner =  cellfun(@isempty,pr2_winners{monkeynum,sessionnum}(1,:));
            model_extension = ['_space',num2str(spacenum),'_within_eval'];

            scatter(...
                avg_pR2.(strcat(model_pairs{1,1},model_extension))(no_winner),...
                avg_pR2.(strcat(model_pairs{1,2},model_extension))(no_winner),...
                [],getColorFromList(1,spacenum+1));
            legend_data(end+1) = scatter(...
                avg_pR2.(strcat(model_pairs{1,1},model_extension))(~no_winner),...
                avg_pR2.(strcat(model_pairs{1,2},model_extension))(~no_winner),...
                [],getColorFromList(1,spacenum+1),'filled');
                
                
            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
            axis square
            if monkeynum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel(sprintf('%s pR2',getModelTitles(model_pairs{1,1})))
            ylabel(sprintf('%s pR2',getModelTitles(model_pairs{1,2})))
        end
        legend(legend_data,spacenames)
    end
    suptitle('Within Pseudo-R^2 pairwise comparisons')
    
    
%% plot, for each model, RT3D vs RT2D predictions on the same plot

    sessionnum = 1;
    legend_data = [];
    figure
    for monkeynum = 1:length(monkey_names)
        for i_mdl = 1:length(models_to_plot) % 1 = RT3D, 2 = RT2D
            % set subplot
            subplot(1,length(monkey_names),monkeynum)
            plot([-1 1],[-1 1],'k--','linewidth',0.5)
            hold on
            plot([0 0],[-1 1],'k-','linewidth',0.5)
            plot([-1 1],[0 0],'k-','linewidth',0.5)
            
            avg_pR2 = neuronAverage(model_eval{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
            % scatter filled circles if there's a winner, empty circles if not
            no_winner =  cellfun(@isempty,pr2_winners{monkeynum,sessionnum}(1,:));
            space_extensions = {'_space1_within_eval','_space2_within_eval'};
            
            legend_data(end+1) = scatter(...
                avg_pR2.(strcat(model_pairs{1,i_mdl},space_extensions{2})),...
                avg_pR2.(strcat(model_pairs{1,i_mdl},space_extensions{1})),...
                [],getColorFromList(1,i_mdl-1),'filled');
            
            %TEMP, plot the link between the two dots in two lines, to see
            %how close the two models' performances are to each other
            %Only plot this in the first iteration
            if i_mdl == 1
                num_neurons = length(avg_pR2.(strcat(model_pairs{1,i_mdl},space_extensions{2})));
                for i = 1:num_neurons
                    model1_x = avg_pR2.(strcat(model_pairs{1,1},space_extensions{2}))(i);
                    model1_y = avg_pR2.(strcat(model_pairs{1,1},space_extensions{1}))(i);
                    model2_x = avg_pR2.(strcat(model_pairs{1,2},space_extensions{2}))(i);
                    model2_y = avg_pR2.(strcat(model_pairs{1,2},space_extensions{1}))(i);
                    plot([model1_x,model2_x],[model1_y,model2_y],'k','LineWidth',2)
                end
            end
                    
                
                
            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-0.1 0.5],'ylim',[-0.1 0.5],'FontSize',16)
            axis square
            if monkeynum ~= 1
                suptitle('Pseudo-R^2 pairwise comparisons')
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel(spacenames(2))
            ylabel(spacenames(1))
        end
        
        models_to_plot_QiweiNames = strings;
        for i = 1:length(models_to_plot)
            if strcmp(string(models_to_plot(i)), 'ext') == 1
                models_to_plot_QiweiNames(i) = 'Hand-Only';
            elseif strcmp(string(models_to_plot(i)), 'handelbow') == 1
                 models_to_plot_QiweiNames(i) = 'Whole-Arm';
            else
                models_to_plot_QiweiNames(i) = 'Unknown Model';
            end
        end
        
        legend(legend_data,models_to_plot_QiweiNames)
        %set(gca,'DefaultTextFontSize',22)
    end
%% (New) Plots between whole-arm and hand-only for RT2D and RT3D
%     
%     task_names = {'Random Target 2D task', 'Random Target 3D task'};
%     % show scatter plot for hand/elbow pR2 within condition vs against condition
%     for modelnum = 1:length(models_to_plot)
%         figure
%         for monkeynum = 1:length(monkey_names)
%             for spacenum = 1:2
%                 % set subplot
%                 subplot(2,length(monkey_names),(spacenum-1)*length(monkey_names)+monkeynum)
% 
%                 % plot lines
%                 plot([-1 1],[-1 1],'k--','linewidth',0.5)
%                 hold on
%                 plot([0 0],[-1 1],'k-','linewidth',0.5)
%                 plot([-1 1],[0 0],'k-','linewidth',0.5)
% 
%                 % plot out each session
%                 for sessionnum = 1:session_ctr(monkeynum)
%                     avg_pR2 = neuronAverage(model_eval{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
%                     scatter(...
%                         avg_pR2.(sprintf('%s_space%d_within_eval',models_to_plot{1},spacenum)),... %Warning, hardcode (1)
%                         avg_pR2.(sprintf('%s_space%d_within_eval',models_to_plot{2},spacenum)),... %Warning, hardcode (2)
%                         [],session_colors(sessionnum,:),'filled')
%                 end
%                 % make axes pretty
%                 set(gca,'box','off','tickdir','out',...
%                     'xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
%                 axis square
%                 if monkeynum ~= 1 || spacenum ~= 1
%                     set(gca,'box','off','tickdir','out',...
%                         'xtick',[],'ytick',[])
%                 end
%                 xlabel(sprintf('Hand-Only model'))
%                 ylabel(sprintf('Whole-Arm model'))
%                 title(task_names(spacenum))
%             end
%         end
%         suptitle('pR^2 Comparison Between Whole-Arm model vs Hand-Only model')
%     end
%    
%     suptitle('Within Pseudo-R^2 pairwise comparisons')
    
%% plot PDs across spaces (actual, as well as predicted by each model)
    sessionnum = 1;
    spacenames = {'RT3D','RT2D'};
    tasknames = {'none','RW'};
    legend_data = [];
    figure
    models_to_plot_temp = models_to_plot;
    models_to_plot_temp{end+1} = 'S1_FR';
    
    for monkeynum = 1:length(monkey_names)
        for mdlnum = 1:numel(models_to_plot_temp)
            % set subplot
            subplot(1,length(monkey_names),monkeynum); hold on;
            plot([-4 4],[-4 4],'k--','linewidth',0.5)
            hold on
            plot([0 0],[-4 4],'k-','linewidth',0.5)
            plot([-4 4],[0 0],'k-','linewidth',0.5)

            keep_mask = strcmpi(model_tuning{monkeynum,sessionnum}.task,tasknames{1});
            space_1_tuning = model_tuning{monkeynum,sessionnum}(keep_mask,:);
            keep_mask = strcmpi(model_tuning{monkeynum,sessionnum}.task,tasknames{2});
            space_2_tuning = model_tuning{monkeynum,sessionnum}(keep_mask,:);
            
            avg_1_tuning = neuronAverage(space_1_tuning,struct('keycols','signalID','do_ci',false));
            avg_2_tuning = neuronAverage(space_2_tuning,struct('keycols','signalID','do_ci',false));

            % scatter filled circles if there's a winner, empty circles if not
            model_extension = ['_velPD'];


            legend_data(end+1) = scatter(...
                avg_1_tuning.(strcat(models_to_plot_temp{mdlnum},model_extension)),...
                avg_2_tuning.(strcat(models_to_plot_temp{mdlnum},model_extension)),...
                [],getColorFromList(1,mdlnum+1),'filled'); hold on;


            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-4 4],'ylim',[-4 4])
            axis square
            if monkeynum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel(sprintf('PD during %s',spacenames{1}))
            ylabel(sprintf('PD during %s',spacenames{2}))
        end
        legend(legend_data,models_to_plot_temp)
    end

%% (Try) to save all the showing figures
FolderName = strcat(base_folderpath, 'analysis-results\');   % Your destination folder
if ~exist(FolderName, 'dir')
   mkdir(FolderName)
end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  %??? What do we have in fighandle, can we have fig title?
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  savefig(fullfile(FolderName, [FigName '.fig']));
  saveas(FigHandle,fullfile(FolderName, [FigName '.png']));
  %saveas(FigHandle,strcat(FolderName, [FigName '.pdf']));
  print(FigHandle, strcat(FolderName,FigName),'-dpdf','-bestfit');
  %savefig(fullfile(FolderName, [FigName '.png']));
end

%% close all figures
close all