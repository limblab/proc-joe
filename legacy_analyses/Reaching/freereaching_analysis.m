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
    clear; clc;
    if ispc
        folderpath = 'D:\Lab\Data\FreeReaching\Han_20210623\neural-data\';
    else
        folderpath = '/data/raeed/project-data/limblab/s1-kinematics';
    end
    
    % load data
    file_info = dir(fullfile(folderpath,'*td*.mat'));

    filenames = horzcat({file_info.name})';
    
    % save directory information (for convenience, since this code takes a
    % while)tw
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
    monkey_names = {'Han'};
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
                % smooth marker data
                td_list{i_td} = smoothSignals(td_list{i_td},struct('signals','dlc_pos'));
                
                % set origin as shoulder position at t=0
                td_list{i_td} = setOriginAsShoulder(td_list{i_td},0); % use fixed position (t=0) or set shoulder as 0 for each data point.                
                % get marker velocity
                td_list{i_td} = getDifferential(td_list{i_td},struct('signals','dlc_pos','alias','dlc_vel'));
                % remove time points where dlc tracking is bad
                dlc_idx = [find((strcmpi(td_list{i_td}.dlc_pos_names,['hand2','_x']))),...
                            find((strcmpi(td_list{i_td}.dlc_pos_names,['elbow1','_x'])))];
%                             find((strcmpi(td_list{i_td}.dlc_pos_names,['hand3','_x']))),...
%                             find((strcmpi(td_list{i_td}.dlc_pos_names,['elbow1','_x']))),...
%                             find((strcmpi(td_list{i_td}.dlc_pos_names,['elbow2','_x']))),...
%                             find((strcmpi(td_list{i_td}.dlc_pos_names,['wrist1','_x']))),...
%                             find((strcmpi(td_list{i_td}.dlc_pos_names,['wrist2','_x']))),...
%                             find((strcmpi(td_list{i_td}.dlc_pos_names,['shoulder1','_x'])))];
%                 bad_points = any(isnan(td_list{i_td}.dlc_pos(:,:)),2) | any(isnan(td_list{i_td}.dlc_vel(:,:)),2);
                bad_points = any(isnan(td_list{i_td}.dlc_pos(:,dlc_idx)),2) | any(isnan(td_list{i_td}.dlc_vel(:,dlc_idx)),2);
                td_names = fieldnames(td_list{i_td});
                for i_name = 1:numel(td_names)
                    if(size(bad_points,1) == size(td_list{i_td}.(td_names{i_name}),1))
                        td_list{i_td}.(td_names{i_name})(bad_points,:) = [];
                    end
                end
                fprintf('Removed %.2f%% percent of trials because of missing markers....',sum(bad_points)/numel(bad_points)*100)
            end
            % remove trials where monkey's arm position is out of the
            % "workspace"
            %td_list{i_td}
            %task_list{i_td}
            [td_list{i_td}, experiment_phase_mask] = getExperimentPhase(td_list{i_td},task_list{i_td});
            fprintf('Removed %.2f%% percent of non_experiment points\n',sum(~experiment_phase_mask)/numel(experiment_phase_mask)*100)
% 
% 
%             % get robot height (z data of a hand marker). This is
%             % meaningless during free reach...
%             markername = 'hand3';
%             dlc_idx = find((strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_z'])));
%             robot_height(end+1) = mean(td_list{i_td}.dlc_pos(:,dlc_idx),'omitnan');

            % overwrite trial-id
            td_list{i_td}.trial_id = (1:1:length(td_list{i_td}.dlc_pos))';
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
        
    
%% loop through files and compute cross correlation between neural signals and kinematics
    for filenum = 1:length(td_all)
        td_list = td_all{filenum};
        task_list = task_list_all{filenum};
        
        % get kinematic data and make plots
        xcorr_data = crossCorrKinNeural(td_list);
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
        
        % Get encoding models, spacenum = 1:2 % 1 = 2D, 2 = 3D
        encoderResults_cell{filenum} = reachingEncoders(td_list,task_list,robot_height_all{filenum},struct(...
            'model_aliases',{included_models},...
            'arrayname',arrayname,...
            'num_tuning_bins',16,...
            'crossval_lookup',[],...
            'get_tuning_curves',true,...
            'num_repeats',20,... % Raeed used 20 repeats.
            'num_folds',5,...
            'subsample','none')); % 'speed','half_data','yaxis','xaxis','zaxis','none'
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
    [model_eval,model_tuning,tuning_corr,shift_vaf,tuned_neurons,model_vars_within] = deal(cell(length(monkey_names),size(session_colors,1)));
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
                encoderResults.crossTuning.(sprintf('glm_%s_model_dlc_vel_handCurve',included_models{modelnum})),...
                encoderResults.crossTuning.(sprintf('glm_%s_model_dlc_vel_handPD',included_models{modelnum})),...
                encoderResults.crossTuning.(sprintf('glm_%s_model_dlc_vel_handPD_zAng',included_models{modelnum})),...
                encoderResults.crossTuning.(sprintf('glm_%s_model_dlc_vel_handModdepth3D',included_models{modelnum})),...
                'VariableNames',strcat(included_models(modelnum),{'_velCurve','_velPD','_velPDzAng','_modDepth'}));
            model_tuning_cell{modelnum}.Properties.VariableDescriptions = {'linear','circular','circular','linear'};
        end
        model_tuning_cell{end} = table(...
            encoderResults.crossTuning.([arrayname,'_FR_dlc_vel_handCurve']),...
            encoderResults.crossTuning.([arrayname,'_FR_dlc_vel_handPD']),...
            encoderResults.crossTuning.([arrayname,'_FR_dlc_vel_handPD_zAng']),...
            encoderResults.crossTuning.([arrayname,'_FR_dlc_vel_handModdepth3D']),...
            'VariableNames',strcat('S1_FR',{'_velCurve','_velPD','_velPDzAng','_modDepth'}));
        model_tuning_cell{end}.Properties.VariableDescriptions = {'linear','circular','circular','linear'};
                
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
    
  
    
%% get predictions across entire file using correct model fold
    td_list = td_all{1};
    for i_td = 1:numel(td_list)
        if(strcmpi(td_list{i_td}.task,'RW'))
            task_idx=1;
        else
            task_idx=2;
        end
        % get predictions

        for mdlnum = 1:numel(included_models)
            % preallocate space
            td_list{i_td}.(encoderResults.params.model_names{mdlnum}) = nan(size(td_list{i_td}.LeftS1_FR,1),size(td_list{i_td}.LeftS1_FR,2),encoderResults.params.num_repeats);
            for i_repeat = 1:encoderResults.params.num_repeats
                for i_fold = 1:encoderResults.params.num_folds
                    % get predictions across whole file
                    temp_td = getModel(td_list{i_td},...
                        encoderResults.glm_info{i_fold +(i_repeat-1)*encoderResults.params.num_folds,mdlnum});

                    % store those that were in the test set. field name =
                    % encoderResults.params.model_names{mdlnum}
                    test_set_idx = encoderResults.params.crossval_lookup.trialID(encoderResults.params.crossval_lookup.spaceNum==task_idx & ...
                        encoderResults.params.crossval_lookup.crossvalID(:,1) == i_repeat & ...
                        encoderResults.params.crossval_lookup.crossvalID(:,2) == i_fold);

                    td_list{i_td}.(encoderResults.params.model_names{mdlnum})(test_set_idx,:,i_repeat) = temp_td.(encoderResults.params.model_names{mdlnum})(test_set_idx,:);
                end
            end
%         
        end
    end

    
    
    
%% %% plot random "trials" and predictions from one of the models - real FR, model predict FR for multiple neurons

    unit_idx_list = [9];
   
    num_rows = length(unit_idx_list) + 1; %last row for speed 
    td_idx = 2; % pick which task to use
    
    window_plot = [-2,2]; % s
    num_trials_plot = 4;
    num_neurons = length(unit_idx_list);
    bin_size = td_list{td_idx}.bin_size;
    % pick random time point, make sure data in window are consecutive
    % based on trial_idx, if yes go on and plot
    
    td_list_smooth = smoothSignals(td_list{td_idx},'LeftS1_FR');

    %randomly choose a time point
    for i_trial = 1:num_trials_plot
        is_consec = 0;
        test_ctr = 0;
        while ~is_consec && test_ctr < 100
            time_idx(i_trial) = (rand()*(size(td_list{td_idx}.dlc_pos,1) - floor(diff(window_plot)*2/bin_size))) + floor(diff(window_plot)/bin_size);
            window_idx(i_trial,:) = time_idx(i_trial)+floor(window_plot./bin_size);
            
            trial_idx(i_trial) = time_idx(i_trial);

            if(max(diff(trial_idx)) >= diff(window_plot))
                is_consec = 1;
            end

            test_ctr = test_ctr + 1;
        end
    end

    window_length = ( window_plot(2) - window_plot(1) )/bin_size;

    f=figure('Position',[10 50 1440 600]); hold on;
    f.Name = [encoderResults.crossEval.monkey{1} ,'_', encoderResults.crossEval.date{1},'_encoderModelExampleFR'];
    % plot FR for unit
    x_data = (0:1:diff(window_idx(1,:)))*bin_size;
    encoderResults = encoderResults_cell{1};
    repeatnum = 1;
    
    for i_neuron = 1:num_rows
        for i_trial = 1:num_trials_plot
            %i_trial
            current_plot_num = (i_neuron - 1) * num_trials_plot + i_trial;
            ax_list(i_trial) = subplot(num_rows,num_trials_plot,current_plot_num); hold on;
            
            %plot the actual firing rate data
            if i_neuron ~= num_rows
                plot(x_data,td_list_smooth.LeftS1_FR(window_idx(i_trial,1):window_idx(i_trial,2),unit_idx_list(i_neuron)),'k','linewidth',2)
                
                % plot the model firing rate data
                plot(x_data,td_list_smooth.glm_ext_model(window_idx(i_trial,1):window_idx(i_trial,2),unit_idx_list(i_neuron),1),'color',getColorFromList(1,2),'linewidth',2)
                plot(x_data,td_list_smooth.glm_handelbow_model(window_idx(i_trial,1):window_idx(i_trial,2),unit_idx_list(i_neuron),1),'color',getColorFromList(1,3),'linewidth',2)
%                 plot(x_data,td_list_smooth.glm_joint_model(window_idx(i_trial,1):window_idx(i_trial,2),unit_idx_list(i_neuron),1),'color',getColorFromList(1,1),'linewidth',2)
                
            end
            
            if i_trial == 1 && i_neuron ~= num_rows%the first column
                ylabel('neuron ' + string(unit_idx_list(i_neuron)));
            end
            
            %Plot speed of hand markers
            if i_neuron == num_rows %if last row
                dlc_section_velocity = td_list_smooth.dlc_vel(window_idx(i_trial,1):window_idx(i_trial,2),19:21);
                dlc_section_speed = sqrt(sum(dlc_section_velocity.^2,2)); % apparently there's a unit conversion here
                plot(x_data,dlc_section_speed,'color','k','linewidth',2);
             
                if i_trial == 1
                    ylabel('speed (mm/s)')
                end
            end

            formatForLee(gcf);
            set(gca,'fontsize',14);
            xlim([0,4])
            %don't plot the x and y ticks if not the first column
            if i_trial ~= 1
                set(gca, 'YTick', [])
            end
            if i_neuron ~= num_rows
                set(gca,'XTick',[])
            end
            
            if(i_trial==1 && i_neuron == 1)
                l=legend('S1','hand-only','whole-arm','joint-ang');
                set(l,'box','off');
            end
        end
    end    
    
    
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
    spacenames = {'RT2D','RT3D'};
    legend_data = [];

    f=figure();
    f.Name = [encoderResults.crossEval.monkey{1} ,'_', encoderResults.crossEval.date{1},'_',model_pairs{1},'_',model_pairs{2},'_encoderModelComparison'];
    for monkeynum = 1:length(monkey_names)
        for spacenum = 1:2 % 1 = 2D, 2 = 3D
            % set subplot
            subplot(1,length(monkey_names),monkeynum)
            plot([-1 1],[-1 1],'k--','linewidth',0.5)
            hold on
            plot([0 0],[-1 1],'k-','linewidth',0.5)
            plot([-1 1],[0 0],'k-','linewidth',0.5)
            
            avg_pR2 = neuronAverage(model_eval{monkeynum,sessionnum},struct('keycols','signalID','do_ci',true));
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
        l=legend(legend_data,spacenames);
%         l=legend(legend_data,'low spd','high spd');
        set(l,'box','off');
    end
    formatForLee(gcf);
    set(gca,'fontsize',14);
    
% plot, for each model, RT3D vs RT2D predictions on the same plot

    sessionnum = 1;
    legend_data = [];
    f=figure;
    f.Name = [encoderResults.crossEval.monkey{1} ,'_', encoderResults.crossEval.date{1},'_',model_pairs{1},'_',model_pairs{2},'_encoderTaskComparison'];
    for monkeynum = 1:length(monkey_names)
        for i_mdl = 1:length(models_to_plot) % 1 = RT2D, 2 = RT3D
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
                avg_pR2.(strcat(model_pairs{1,i_mdl},space_extensions{1})),...
                avg_pR2.(strcat(model_pairs{1,i_mdl},space_extensions{2})),...
                [],getColorFromList(1,i_mdl-1),'filled');
            
            % link corresponding data points.
            if i_mdl == 1
                num_neurons = length(avg_pR2.(strcat(model_pairs{1,i_mdl},space_extensions{2})));
                for i = 1:num_neurons
                    model1_x = avg_pR2.(strcat(model_pairs{1,1},space_extensions{1}))(i);
                    model1_y = avg_pR2.(strcat(model_pairs{1,1},space_extensions{2}))(i);
                    model2_x = avg_pR2.(strcat(model_pairs{1,2},space_extensions{1}))(i);
                    model2_y = avg_pR2.(strcat(model_pairs{1,2},space_extensions{2}))(i);
                    plot([model1_x,model2_x],[model1_y,model2_y],'k','LineWidth',2)
                end
            end
    
            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-0.1 0.5],'ylim',[-0.1 0.5])
            axis square
            xlabel(spacenames(1))
            ylabel(spacenames(2))
%             xlabel('RT3D low spd pR2');
%             ylabel('RT3D high spd pR2');
        end
        
        models_to_plot_names = strings;
        for i = 1:length(models_to_plot)
            models_to_plot_names(i) = getModelTitles(model_pairs{1,i});
        end
        
        l=legend(legend_data,models_to_plot_names);
        set(l,'box','off');
        set(gca,'FontSize',14)
        formatForLee(gcf);

    end
    
    
%% histogram of whole-arm pR2 - hand-elbow pR2 for each space, pair data

    sessionnum = 1;
    spacenames = {'RT2D','RT3D'};
    legend_data = [];
    data = [];
    
    pr2_diff_winners = cell(length(monkey_names),size(session_colors,1));
    for monkeynum = 1:length(monkey_names)
        for sessionnum = 1:session_ctr(monkeynum)
            data = [];
            for spacenum = 1:2
                model_extension = ['_space',num2str(spacenum),'_within_eval'];
                % store difference metric in model_eval table
                data(:,spacenum) = model_eval{monkeynum,sessionnum}.(strcat(model_pairs{1,2},model_extension)) - model_eval{monkeynum,sessionnum}.(strcat(model_pairs{1,1},model_extension));
                model_eval{monkeynum,sessionnum}.(['model_comp_diff',model_extension]) = data(:,spacenum);
            end
            
            % get winners
            [pr2_diff_winners{monkeynum,sessionnum}] = compareEncoderMetrics(...
                    model_eval{monkeynum,sessionnum},struct(...
                        'bonferroni_correction',1,...
                        'models',{models_to_plot},...
                        'model_pairs',{{'model_comp_diff_space1','model_comp_diff_space2'}},...
                        'postfix','_within_eval'));
                    
            no_winner =  cellfun(@isempty,pr2_diff_winners{monkeynum,sessionnum}(1,:));
        end
    end
    
    data = [];
    for monkeynum = 1:length(monkey_names)
        for spacenum = 1:2 % 1 = 2D, 2 = 3D
            model_extension = ['_space',num2str(spacenum),'_within_eval'];
            avg_pR2 = neuronAverage(model_eval{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
            % scatter filled circles if there's a winner, empty circles if not
%             no_winner =  cellfun(@isempty,pr2_winners{monkeynum,sessionnum}(1,:));
            

            data(:,spacenum) = avg_pR2.(strcat(model_pairs{1,2},model_extension)) - avg_pR2.(strcat(model_pairs{1,1},model_extension));
        end
    end
    f=figure();
    f.Name = [encoderResults.crossEval.monkey{1} ,'_', encoderResults.crossEval.date{1},'_',model_pairs{1},'_',model_pairs{2},'_pr2performancecomparison'];
    f.Position = [680 485 1058 493];
    subplot(1,2,1); hold on;
    
    plot(data(~no_winner,1),data(~no_winner,2),'.','color',getColorFromList(1,1),'markersize',20)
    plot(data(no_winner,1),data(no_winner,2),'o','color',getColorFromList(1,1),'markersize',6)
    
    xlabel('RT2D whole-arm pR^2 minus hand-only pR^2');
    ylabel('RT3D whole-arm pR^2 minus hand-only pR^2');
    plot([-0.1,0.5],[-0.1,0.5],'k--');
%     xlim([-0.02,0.1]);
%     ylim([-0.02,0.1]);
    formatForLee(gcf);
    set(gca,'fontsize',14);
    
    subplot(1,2,2)
    histogram(data(:,2)-data(:,1),-0.1:0.005:0.1)
    formatForLee(gcf);
    set(gca,'fontsize',14);
    xlabel('Whole-arm pR^2 - Hand-only pR^2 across speeds')
    ylabel('Number of neurons');
    [p] = ranksum(data(:,2),data(:,1),'tail','right')
    

%% get predictions during high velocity periods. Compare predicted FR to actual for each model
% for each task
mse = zeros(size(td_list{1}.LeftS1_spikes,2),2,numel(td_list));
for i_task = 1:numel(td_list)
    if(strcmpi(td_list{i_td}.task,'RW'))
        space_idx = 1; % this corresponds to the space idx
    else
        space_idx=0;
    end
    
    % get peaks in speed
    dlc_spd = sqrt(sum(td_list{i_task}.dlc_vel.^2,2));
    [pks, locs] = findpeaks(dlc_spd,'MinPeakHeight',prctile(dlc_spd,85),'MinPeakDistance',10);
%     figure();
%     plot(dlc_spd)
%     hold on
%     plot(locs, prctile(dlc_spd,90)*ones(size(locs)),'.','markersize',12)
    
    % get FR, predicted FR for each neuron at peaks in speed
    FR_real = td_list{i_task}.LeftS1_spikes(locs,:);
    FR_ext = mean(td_list{i_task}.glm_ext_model(locs,:,:),3,'omitnan');
    FR_whole = mean(td_list{i_task}.glm_handelbow_model(locs,:,:),3,'omitnan'); % average over repeats
    
    mse(:,1,i_task) = mean((FR_ext-FR_real).^2,1,'omitnan');
    mse(:,2,i_task) = mean((FR_whole-FR_real).^2,1,'omitnan');
%     
%     figure();
%     for i_unit = 1:size(FR_real,2)
%         subplot(5,5,i_unit)
%         plot(FR_real(:,i_unit), FR_ext(:,i_unit),'.')
%         hold on
%         plot(FR_real(:,i_unit), FR_whole(:,i_unit),'o','markersize',4)
%     end
end

mse_model_diff = squeeze(mse(:,2,:) - mse(:,1,:));
    
%% get FR during reaches in the PD for each neuron
    td_list = td_all{1};
    task_list = task_list_all{1};
    cutoff = 95;
    
    tasknames = {'RW','none'};
    keep_mask = strcmpi(model_tuning{1,1}.task,tasknames{1});
    space_1_tuning = model_tuning{1,1}(keep_mask,:);
    keep_mask = strcmpi(model_tuning{1,1}.task,tasknames{2});
    space_2_tuning = model_tuning{1,1}(keep_mask,:);
    avg_1_tuning = neuronAverage(space_1_tuning,struct('keycols','signalID','do_ci',true,'do_nanmean',true));
    avg_2_tuning = neuronAverage(space_2_tuning,struct('keycols','signalID','do_ci',true,'do_nanmean', true));

    td_2d = td_list{find(strcmpi(task_list,'RT'))};
    td_3d = td_list{find(strcmpi(task_list,'RT3D'))};
    is_3d = nan;
    FR_high_spd = [];
    for i_unit = 1:length(td_list{1}.LeftS1_unit_guide)
        avg_1_idx = find(all(avg_1_tuning.signalID == td_list{1}.LeftS1_unit_guide(i_unit,:),2));
        velPD_2D = avg_1_tuning.S1_FR_velPD(avg_1_idx);
        velPD_zAng_2D = avg_1_tuning.S1_FR_velPDzAng(avg_1_idx);
        
        avg_2_idx = find(all(avg_2_tuning.signalID == td_list{1}.LeftS1_unit_guide(i_unit,:),2));
        velPD_3D = avg_2_tuning.S1_FR_velPD(avg_2_idx);
        velPD_zAng_3D = avg_2_tuning.S1_FR_velPDzAng(avg_2_idx);

        PD_2D_vec = convertToVec(velPD_2D, velPD_zAng_2D);
        PD_2D_vec = PD_2D_vec(1:2);
        PD_3D_vec = convertToVec(velPD_3D, velPD_zAng_3D);
        
        % for each trial data, get high speed movements in direction of PD
        % then get firing rates during those movements
        for i_td = 1:2
            switch i_td
                case 1
                    td_use = td_2d;
                    vec_use = PD_2D_vec;
                    is_3d = 0;
                case 2
                    td_use = td_3d;
                    vec_use = PD_3D_vec;
                    is_3d = 1;
            end
            dlc_idx = [find((strcmpi(td_use.dlc_pos_names,['hand2','_x']))),...
                            find((strcmpi(td_use.dlc_pos_names,['hand2','_y']))),...
                           find((strcmpi(td_use.dlc_pos_names,['hand2','_z'])))];
                        
            pd_len = numel(vec_use);
            spd_data = sqrt(sum(td_use.dlc_vel(:,dlc_idx(1:pd_len)).^2,2));
            
            cos_sim = td_use.dlc_vel(:,dlc_idx(1:pd_len))*vec_use'./sqrt(sum(vec_use.^2,2))./spd_data;
            
            in_dir = cos_sim > cos(pi/4);
            FR_data = td_use.LeftS1_FR(in_dir,i_unit);
            FR_data(isnan(FR_data)) = [];
            
            FR_high_spd(i_unit,i_td) = prctile(FR_data,cutoff);
            
        end
    end
    
    figure();
    bin_edges = -95:10:95;
    histogram(FR_high_spd(:,2)-FR_high_spd(:,1),bin_edges);

%% save variables
    % 3D PD
    tasknames = {'RW','none'};
    keep_mask = strcmpi(model_tuning{1,1}.task,tasknames{1});
    space_1_tuning = model_tuning{1,1}(keep_mask,:);
    keep_mask = strcmpi(model_tuning{1,1}.task,tasknames{2});
    space_2_tuning = model_tuning{1,1}(keep_mask,:);
    avg_1_tuning = neuronAverage(space_1_tuning,struct('keycols','signalID','do_ci',true,'do_nanmean',true));
    avg_2_tuning = neuronAverage(space_2_tuning,struct('keycols','signalID','do_ci',true,'do_nanmean', true));


    tuning_data = [];
    tuning_data.velPD_2D = avg_1_tuning.S1_FR_velPD;
    tuning_data.velPD_2D_CILo = avg_1_tuning.S1_FR_velPDCILo;
    tuning_data.velPD_2D_CIHi = avg_1_tuning.S1_FR_velPDCIHi;
    
    tuning_data.velPD_zAng_2D = avg_1_tuning.S1_FR_velPDzAng;
    tuning_data.velPD_zAng_2D_CILo = avg_1_tuning.S1_FR_velPDzAngCILo;
    tuning_data.velPD_zAng_2D_CIHi = avg_1_tuning.S1_FR_velPDzAngCIHi;
    
    tuning_data.velPD_3D = avg_2_tuning.S1_FR_velPD;
    tuning_data.velPD_3D_CILo = avg_2_tuning.S1_FR_velPDCILo;
    tuning_data.velPD_3D_CIHi = avg_2_tuning.S1_FR_velPDCIHi;
    
    tuning_data.velPD_zAng_3D = avg_2_tuning.S1_FR_velPDzAng;
    tuning_data.velPD_zAng_3D_CILo = avg_2_tuning.S1_FR_velPDzAngCILo;
    tuning_data.velPD_zAng_3D_CIHi = avg_2_tuning.S1_FR_velPDzAngCIHi;

      
%% save all the showing figures

% fpath = folderpath; % folderpath is defined at the beginning, can overwrite to a different path

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  %
  if(strcmpi(FigHandle.Name, '')==1)
    FigHandle.Name = num2str(get(FigHandle, 'Number'));
  end
  saveFiguresLIB(FigHandle, fpath, FigHandle.Name);

end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% get pR2 vs. speed for both models
    td_list = td_all{1};
    task_idx = find(strcmpi(task_list,'RT3D'));
    
    % get speed
    markername = 'hand2';
    dlc_idx = [find((strcmpi(td_list{task_idx}.dlc_pos_names,[markername,'_x']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,[markername,'_y']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,[markername,'_z'])))];
    
    spd_data = sqrt(sum(td_list{task_idx}.dlc_vel(:,dlc_idx).^2,2));

    % get predictions
    
    for mdlnum = 1:numel(included_models)
        % preallocate space
        td_list{task_idx}.(encoderResults.params.model_names{mdlnum}) = nan(size(td_list{task_idx}.LeftS1_FR,1),size(td_list{task_idx}.LeftS1_FR,2),encoderResults.params.num_repeats);
        for i_repeat = 1:encoderResults.params.num_repeats
            for i_fold = 1:encoderResults.params.num_folds
                % get predictions across whole file
                temp_td = getModel(td_list{task_idx},...
                    encoderResults.glm_info{i_fold +(i_repeat-1)*encoderResults.params.num_folds,mdlnum});
                
                % store those that were in the test set. field name =
                % encoderResults.params.model_names{mdlnum}
                test_set_idx = encoderResults.params.crossval_lookup.trialID(encoderResults.params.crossval_lookup.spaceNum==2 & ...
                    encoderResults.params.crossval_lookup.crossvalID(:,1) == i_repeat & ...
                    encoderResults.params.crossval_lookup.crossvalID(:,2) == i_fold);
                
                td_list{task_idx}.(encoderResults.params.model_names{mdlnum})(test_set_idx,:,i_repeat) = temp_td.(encoderResults.params.model_names{mdlnum})(test_set_idx,:);
            end
        end
%         
    end
    
    %% compute prediction error vs speed
    mean_pred_error = zeros(size(td_list{task_idx}.LeftS1_FR,1),size(td_list{task_idx}.LeftS1_FR,2),2);
    for mdlnum = 1:numel(included_models)  
        mean_pred_error(:,:,mdlnum) = mean(abs(td_list{task_idx}.LeftS1_FR - td_list{task_idx}.(encoderResults.params.model_names{mdlnum})),3,'omitnan');
    end
    
    figure()
    idx = 6;
    plot(spd_data,mean_pred_error(:,idx,2) - mean_pred_error(:,idx,1),'.')
    
    
%% get pR2 pairwise and histogram diff for high correlation conditions   
    

    % evaluate model on segments of data. Compare performance vs.
    % correlation
    seg_len = 20; % 50 ms bins
    seg_start = 1:seg_len:(size(td_list{task_idx}.LeftS1_FR,1)-seg_len);
    
    seg_end = seg_start + seg_len-1;
    temp_td = td_list{task_idx};
    n_repeats = encoderResults.params.num_repeats;
    
    dlc_idx = [find((strcmpi(td_list{task_idx}.dlc_pos_names,['hand2_x']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,['hand2_y']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,['elbow1_x']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,['elbow1_y'])))];
    
    eval_ext = zeros(numel(seg_start),size(td_list{task_idx}.LeftS1_FR,2),n_repeats);
    eval_handelbow = zeros(size(eval_ext));
    corr_data = zeros(numel(seg_start),2); % hand-elbow x, hand-elbow y correlation
    for i_seg = 1:numel(seg_start)    
        for i_rep = 1:n_repeats
            % get pr2 vals for this segment
            temp_td.LeftS1_FR = td_list{task_idx}.LeftS1_FR(seg_start(i_seg):seg_end(i_seg),:);
            temp_td.glm_ext_model = td_list{task_idx}.glm_ext_model(seg_start(i_seg):seg_end(i_seg),:,i_rep);
            temp_td.glm_handelbow_model = td_list{task_idx}.glm_handelbow_model(seg_start(i_seg):seg_end(i_seg),:,i_rep);
            eval_ext(i_seg,:,i_rep) = evalModel(temp_td,struct('out_signals','LeftS1_FR','model_name','ext_model','model_type','glm','eval_metric','pr2'));
            eval_handelbow(i_seg,:,i_rep) = evalModel(temp_td,struct('out_signals','LeftS1_FR','model_name','handelbow_model','model_type','glm','eval_metric','pr2'));
        end
        % get correlation values for this segment
        temp_corr = corr(td_list{task_idx}.dlc_vel(seg_start(i_seg):seg_end(i_seg),dlc_idx));
        corr_data(i_seg,:) = [temp_corr(1,3),temp_corr(2,4)];
    end
        
    mean_eval_ext = mean(eval_ext,3);
    mean_eval_handeblow = mean(eval_handelbow,3);
       
    histogram(mean_eval_handeblow - mean_eval_ext)
    
    
    
    
    
%% plot vel and pos magnitude against each other for all neurons
    % need to account for magnitude of data when doing this
    % only use extrinsic (hand-only) model
    
    % get for each repeat and fold (to measure spread)
    
    model_use = 'ext_model';
    
    % get model idx
    for i=1:size(encoderResults.glm_info_within,3)
        if(strcmpi(encoderResults.glm_info_within{1,1,i}.model_name,model_use)==1)
            glm_idx = i;
        end
    end
    
    % normalize by magnitude of corresponding data
    spacename = 'RT3D'; % space name sets range of velocities and positions
    td_idx = find(strcmpi(task_list,spacename),1,'first');
    pos_data = td_list{td_idx}.dlc_pos(:,encoderResults.glm_info_within{1, 2,glm_idx}.in_signals{1,2});
    vel_data = td_list{td_idx}.dlc_vel(:,encoderResults.glm_info_within{1, 2,glm_idx}.in_signals{2,2});
    
    pos_range = prctile(pos_data,95) - prctile(pos_data,5);
    vel_range = prctile(vel_data,95) - prctile(vel_data,5);
        
    pos_norm = sqrt(sum(pos_range.^2));
    vel_norm = sqrt(sum(vel_range.^2));
    
    % get normalized glm_param vals
    b_list = zeros(size(encoderResults.glm_info_within,1),6,size(encoderResults.glm_info_within{1,1,1}.b,2));% hand pos (x,y,z) vel (x,y,z)
    proj_list = zeros(size(b_list,1),2,size(b_list,3));
    for i_rep = 1:size(encoderResults.glm_info_within,1)
        b_list(i_rep,:,:) = encoderResults.glm_info_within{i_rep, glm_idx}.b(2:end,:);
        
        %glm_info_within(,1,) is RT2d; 2 = RT3d
        proj_list(i_rep,1,:) = sqrt(sum((encoderResults.glm_info_within{i_rep, 2,glm_idx}.b(2:4,:)).^2,1)); % pos proj
        proj_list(i_rep,2,:) = sqrt(sum((encoderResults.glm_info_within{i_rep, 2,glm_idx}.b(5:end,:)).^2,1)); % vel proj
    end
    
    mean_proj = squeeze(mean(proj_list,1)).*[pos_norm; vel_norm];
    
    
    
    
    max_lim = max(max(mean_proj))*1.1;
    figure();
    plot(mean_proj(1,:),mean_proj(2,:),'.','color',getColorFromList(1,0),'markersize',12)
    hold on
    plot([0,max_lim],[0,max_lim],'k--','linewidth',2);
    
    xlabel('Pos weight');
    ylabel('Vel weight');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    xlim([0,max_lim]);
    ylim([0,max_lim]);
    
    
%% plot mod depth across tasks
    figure();
    
    sessionnum=1;
    tasknames = {'RW','none'}; % RT2d = 1, RT3d = 2
    for monkeynum = 1:length(monkey_names)
        for spacenum = 1:2 % 1 = 2D, 2 = 3D
            % set subplot
            subplot(1,length(monkey_names),monkeynum)
            plot([-1 1],[-1 1],'k--','linewidth',0.5)
            hold on
            plot([0 0],[-1 1],'k-','linewidth',0.5)
            plot([-1 1],[0 0],'k-','linewidth',0.5)
            
            keep_mask = strcmpi(model_tuning{monkeynum,sessionnum}.task,tasknames{1});
            space_1_tuning = model_tuning{monkeynum,sessionnum}(keep_mask,:);
            keep_mask = strcmpi(model_tuning{monkeynum,sessionnum}.task,tasknames{2});
            space_2_tuning = model_tuning{monkeynum,sessionnum}(keep_mask,:);
            
            avg_1_tuning = neuronAverage(space_1_tuning,struct('keycols','signalID','do_nanmean',true));
            avg_2_tuning = neuronAverage(space_2_tuning,struct('keycols','signalID','do_nanmean', true));
            % scatter filled circles if there's a winner, empty circles if not
           
            plot(avg_1_tuning.S1_FR_modDepth,avg_2_tuning.S1_FR_modDepth,'.')
                
                
            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
            axis square
            if monkeynum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel('RT2D Mod Depth');
            ylabel('RT3D Mod Depth');
        end
    end
    formatForLee(gcf);
    set(gca,'fontsize',14);

    
    
%% plot PDs across spaces (bootstraps based on actual FR)
    sessionnum = 1;
    spacenames = {'RT2D','RT3D'};
    tasknames = {'RW','none'};
    legend_data = [];
    f=figure();
    f.Name = [encoderResults.crossEval.monkey{1} ,'_', encoderResults.crossEval.date{1},'_PDTaskComparison'];
    models_to_plot_temp = {'S1_FR'};
    max_lim = 2*pi;
    for monkeynum = 1:length(monkey_names)
        for mdlnum = 1:numel(models_to_plot_temp)
            % set subplot
            subplot(1,length(monkey_names),monkeynum); hold on;
            plot([-max_lim max_lim],[-max_lim max_lim],'k--','linewidth',0.5)
            plot([-max_lim max_lim],pi/2+[-max_lim max_lim],'k:','linewidth',0.5)
            plot([-max_lim max_lim],-pi/2+[-max_lim max_lim],'k:','linewidth',0.5)
            hold on
            plot([0 0],[-max_lim max_lim],'k-','linewidth',0.5)
            plot([-max_lim max_lim],[0 0],'k-','linewidth',0.5)

            keep_mask = strcmpi(model_tuning{monkeynum,sessionnum}.task,tasknames{1});
            space_1_tuning = model_tuning{monkeynum,sessionnum}(keep_mask,:);
            keep_mask = strcmpi(model_tuning{monkeynum,sessionnum}.task,tasknames{2});
            space_2_tuning = model_tuning{monkeynum,sessionnum}(keep_mask,:);
            
            avg_1_tuning = neuronAverage(space_1_tuning,struct('keycols','signalID','do_ci',true,'do_nanmean',true));
            avg_2_tuning = neuronAverage(space_2_tuning,struct('keycols','signalID','do_ci',true,'do_nanmean', true));

            % scatter filled circles if there's a winner, empty circles if not
            model_extension = ['_velPD'];

            % deal with circular data
            avg_1_pd = avg_1_tuning.(strcat(models_to_plot_temp{mdlnum},model_extension));
            avg_1_ci = [avg_1_tuning.(strcat(models_to_plot_temp{mdlnum},model_extension,'CILo')), avg_1_tuning.(strcat(models_to_plot_temp{mdlnum},model_extension,'CIHi'))];
            avg_2_pd = avg_2_tuning.(strcat(models_to_plot_temp{mdlnum},model_extension));
            avg_2_ci = [avg_2_tuning.(strcat(models_to_plot_temp{mdlnum},model_extension,'CILo')), avg_2_tuning.(strcat(models_to_plot_temp{mdlnum},model_extension,'CIHi'))];
            
            avg_2_ci = avg_1_pd + angleDiff(avg_1_pd,avg_2_ci,1,1); % use radians, do preserve sign
            avg_2_pd = avg_1_pd + angleDiff(avg_1_pd,avg_2_pd,1,1); % use radians, do preserve sign
            
            legend_data(end+1) = errorbar(avg_1_pd,avg_2_pd,...
                abs(angleDiff(avg_2_ci(:,1),avg_2_pd,1,1)),abs(angleDiff(avg_2_ci(:,2),avg_2_pd,1,1)),... % yneg, ypos
                abs(angleDiff(avg_1_ci(:,1),avg_2_pd,1,1)),abs(angleDiff(avg_1_ci(:,2),avg_2_pd,1,1)),... % xneg, xpos
                'color',getColorFromList(1,mdlnum),'marker','.','linestyle','none','markersize',20,'linewidth',1.5); hold on;


            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-max_lim max_lim],'ylim',[-max_lim max_lim],'fontsize',14)
            axis square
            xlabel(sprintf('Planar PD during %s',spacenames{1}))
            ylabel(sprintf('Planar PD during %s',spacenames{2}))
        end
    end

%% plot z_ang of PD for 3D task as a histogram
    sessionnum = 1;
    taskname = 'none';
    models_to_plot_temp = {'S1_FR'};
    f=figure();
    f.Name = [encoderResults.crossEval.monkey{1} ,'_', encoderResults.crossEval.date{1},'_PDZAngDistribution'];
    for monkeynum = 1:length(monkey_names)
        for mdlnum = 1:numel(models_to_plot_temp)
            % set subplot
            subplot(1,length(monkey_names),monkeynum); hold on;
      
            keep_mask = strcmpi(model_tuning{monkeynum,sessionnum}.task,taskname);
            space_2_tuning = model_tuning{monkeynum,sessionnum}(keep_mask,:);
                        avg_2_tuning = neuronAverage(space_2_tuning,struct('keycols','signalID','do_ci',true));

            % scatter filled circles if there's a winner, empty circles if not
            model_extension = ['_velPDzAng'];

            avg_2_pd = avg_2_tuning.(strcat(models_to_plot_temp{mdlnum},model_extension));
            avg_2_ci = [avg_2_tuning.(strcat(models_to_plot_temp{mdlnum},model_extension,'CILo')), avg_2_tuning.(strcat(models_to_plot_temp{mdlnum},model_extension,'CIHi'))];
            
            histogram(avg_2_pd,20);
            xlim([-pi/2,pi/2]);
            xlabel('Z-ang (radians)');
            ylabel('Number of neurons');
            formatForLee(gcf);
            set(gca,'fontsize',14);
        end
    end
    
%% plot PD shifts across spaces. Compare model to actual shift
    sessionnum = 1;
    spacenames = {'RT2D','RT3D'};
    tasknames = {'RW','none'};
    legend_data = [];
    f=figure();
    f.Name = [encoderResults.crossEval.monkey{1} ,'_', encoderResults.crossEval.date{1},'_PDModelShift'];
    models_to_compare = {'ext','handelbow'};
    base_model = 'S1_FR';
    % scatter filled circles if there's a winner, empty circles if not
    model_extension = ['_velPD'];
    max_lim = 2*pi;
    for monkeynum = 1:length(monkey_names)
        % get base model shift
        keep_mask = strcmpi(model_tuning{monkeynum,sessionnum}.task,tasknames{1});
        space_1_tuning = model_tuning{monkeynum,sessionnum}(keep_mask,:);
        keep_mask = strcmpi(model_tuning{monkeynum,sessionnum}.task,tasknames{2});
        space_2_tuning = model_tuning{monkeynum,sessionnum}(keep_mask,:);

        avg_1_tuning = neuronAverage(space_1_tuning,struct('keycols','signalID','do_ci',false));
        avg_2_tuning = neuronAverage(space_2_tuning,struct('keycols','signalID','do_ci',false));
        
        avg_1_pd = avg_1_tuning.(strcat(base_model,model_extension));
        avg_2_pd = avg_2_tuning.(strcat(base_model,model_extension));
        avg_base_tuning_shift = angleDiff(avg_1_pd,avg_2_pd,1,1); % use radians, preserve sign.
        
        for mdlnum = 1:numel(models_to_compare)
            % set subplot
            subplot(1,length(monkey_names),monkeynum); hold on;
            plot([-max_lim max_lim],[-max_lim max_lim],'k--','linewidth',0.5)
            plot([-max_lim max_lim],pi/2+[-max_lim max_lim],'k:','linewidth',0.5)
            plot([-max_lim max_lim],-pi/2+[-max_lim max_lim],'k:','linewidth',0.5)
            hold on
            plot([0 0],[-max_lim max_lim],'k-','linewidth',0.5)
            plot([-max_lim max_lim],[0 0],'k-','linewidth',0.5)

            avg_1_pd = avg_1_tuning.(strcat(models_to_compare{mdlnum},model_extension));
            avg_2_pd = avg_2_tuning.(strcat(models_to_compare{mdlnum},model_extension));
            
            avg_model_tuning_shift = angleDiff(avg_1_pd,avg_2_pd,1,1); % use radians, preserve sign.
    
            model_to_plot = avg_base_tuning_shift + angleDiff(avg_base_tuning_shift,avg_model_tuning_shift, 1,1);
        
            legend_data(end+1) = scatter(...
                avg_base_tuning_shift,...
                model_to_plot,...
                [],getColorFromList(1,mdlnum+1),'filled'); hold on;


            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-max_lim max_lim],'ylim',[-max_lim max_lim],'fontsize',14)
            axis square
            if monkeynum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel('Actual Planar PD shift')
            ylabel('Model Planar PD shift')
        end
        l=legend(legend_data,models_to_compare);
        set(l,'box','off');
        formatForLee(gcf);
    end
    
    
    
%% plot z angle PD actual vs. predicted for RT3D
    sessionnum = 1;
    spacenames = {'RT2D','RT3D'};
    tasknames = {'RW','none'};
    legend_data = [];
    f=figure();
    f.Name = [encoderResults.crossEval.monkey{1} ,'_', encoderResults.crossEval.date{1},'_PDZAngModel'];
    models_to_compare = {'ext','handelbow'};
    base_model = 'S1_FR';
    % scatter filled circles if there's a winner, empty circles if not
    model_extension = ['_velPDzAng'];
    max_lim = pi/2;
    for monkeynum = 1:length(monkey_names)
        % get base model shift

        keep_mask = strcmpi(model_tuning{monkeynum,sessionnum}.task,tasknames{2});
        space_2_tuning = model_tuning{monkeynum,sessionnum}(keep_mask,:);

        avg_2_tuning = neuronAverage(space_2_tuning,struct('keycols','signalID','do_ci',false));
        
        avg_2_pd = avg_2_tuning.(strcat(base_model,model_extension));
        
        for mdlnum = 1:numel(models_to_compare)
            % set subplot
            subplot(1,length(monkey_names),monkeynum); hold on;
            plot([-max_lim max_lim],[-max_lim max_lim],'k--','linewidth',0.5)
            hold on
            plot([0 0],[-max_lim max_lim],'k-','linewidth',0.5)
            plot([-max_lim max_lim],[0 0],'k-','linewidth',0.5)

            avg_2_pd_model = avg_2_tuning.(strcat(models_to_compare{mdlnum},model_extension));
                
            model_to_plot = avg_2_pd + angleDiff(avg_2_pd,avg_2_pd_model, 1,1);
        
            legend_data(end+1) = scatter(...
                avg_2_pd,...
                model_to_plot,...
                [],getColorFromList(1,mdlnum+1),'filled'); hold on;


            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-max_lim max_lim],'ylim',[-max_lim max_lim],'fontsize',14)
            axis square
            if monkeynum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel('Actual Z-ang of PD')
            ylabel('Model Z-ang of PD')
        end
        l=legend(legend_data,models_to_compare);
        set(l,'box','off');
        formatForLee(gcf);
    end
    
    
%% compute correlation matrices between neurons for 2D and 3D task

    spacenames = {'RT','RT3D'};
    tree_perm_use = []; max_corr = -100;
    for i_space = 1:numel(spacenames)
        % get td from td_list based on task_list
        space_idx = find(strcmpi(spacenames{i_space}, task_list)==1,1,'first');
        td_use = td_list{space_idx};
        
        % get correlations across neurons
        corr_mat = corr(td_use.LeftS1_FR);
        
        % cluster correlations  if correct task
        if(strcmpi(spacenames{i_space},'RT')==1)
            T = linkage(corr_mat);
            f_to_close = figure();
            [~,~,tree_perm_use] = dendrogram(T);
            close(f_to_close);
        end
        % get max corr for color limits later
        for i = 1:size(corr_mat,1) corr_mat(i,i) = nan; end % set diagonal as nan
        max_corr = max(max_corr,max(max(corr_mat)));
    end
        
    for i_space = 1:numel(spacenames)
        % get td from td_list based on task_list
        space_idx = find(strcmpi(spacenames{i_space}, task_list)==1,1,'first');
        td_use = td_list{space_idx};
        % plot correlations
        f=figure();
        f.Name = [encoderResults.crossEval.monkey{1} ,'_', encoderResults.crossEval.date{1},'_',task_list{space_idx},'_NeuronCorr'];
        corr_mat_use = corr(td_use.LeftS1_FR(:,tree_perm_use));
        for i = 1:size(corr_mat_use,1) corr_mat_use(i,i) = nan; end % set diagonal as nan
        imagesc(corr_mat_use);
        b=colorbar;
        caxis(max_corr*[-1,1]);
        cmap = colorcet('D1'); colormap(cmap);
    end




%% split td based on speed. 1 =  low speed (called RT2D in figs), 2 = high speed (called RT3D in figs)
    
    task_idx = 1;
    td_list = {};
    td_use = td_all{1}{task_idx};
    
    task_list_all{1} = {'RT';'RT3D'};
    
    markername = 'hand2';
    dlc_idx = [find((strcmpi(td_use.dlc_pos_names,[markername,'_x']))),...
                find((strcmpi(td_use.dlc_pos_names,[markername,'_y']))),...
                find((strcmpi(td_use.dlc_pos_names,[markername,'_z'])))];
            
    marker_vel = td_use.dlc_vel(:,dlc_idx);
    marker_spd = sqrt(sum(marker_vel.^2,2));
    
    marker_spd_cutoff = prctile(marker_spd,50);
    low_spd_mask = marker_spd < marker_spd_cutoff;
        
    td_list{1} = td_use;
    td_list{2} = td_use;
    
    td_fields = fieldnames(td_use);
    pos_length = length(td_use.dlc_pos);

    for i_field = 1:numel(td_fields)
        if(size(td_use.(td_fields{i_field}),1) == pos_length)
            td_list{1}.(td_fields{i_field}) = td_list{1}.(td_fields{i_field})(low_spd_mask==1,:);
            td_list{2}.(td_fields{i_field}) = td_list{2}.(td_fields{i_field})(low_spd_mask==0,:);
        end
    end

    td_all{1} = td_list;
    %%
        % simple -- compare predictions for movements in x direction vs
    % movements in y direction
    
    move_dir = atan2(td_list{task_idx}.dlc_vel(:,dlc_idx(2)), td_list{task_idx}.dlc_vel(:,dlc_idx(1)));

    in_y_dir = angleDiff(move_dir, pi/2,1,0) < pi/4 & spd_data > prctile(spd_data,40);
    in_x_dir = angleDiff(move_dir,0,1,0) < pi/4 & spd_data > prctile(spd_data,40);
    
    % get x and y dir predictions and pr2 score
    
    n_repeats = size(td_list{task_idx}.glm_ext_vel_model,3);
    eval_ext_y = zeros(size(td_list{task_idx}.LeftS1_FR,2),n_repeats);
    eval_handelbow_y = zeros(size(eval_ext_y));
    eval_ext_x = zeros(size(eval_ext_y));
    eval_handelbow_x = zeros(size(eval_ext_y));
    temp_td = td_list{task_idx};
    

    for i_rep = 1:n_repeats
        % y dir
        temp_td.LeftS1_FR = td_list{task_idx}.LeftS1_FR(in_y_dir,:);
        temp_td.glm_ext_model = td_list{task_idx}.glm_ext_vel_model(in_y_dir,:,i_rep);
        temp_td.glm_handelbow_model = td_list{task_idx}.glm_handelbow_vel_model(in_y_dir,:,i_rep);
        eval_ext_y(:,i_rep) = evalModel(temp_td,struct('out_signals','LeftS1_FR','model_name','ext_model','model_type','glm','eval_metric','pr2'));
        eval_handelbow_y(:,i_rep) = evalModel(temp_td,struct('out_signals','LeftS1_FR','model_name','handelbow_model','model_type','glm','eval_metric','pr2'));
        % x dir
        temp_td.LeftS1_FR = td_list{task_idx}.LeftS1_FR(in_x_dir,:);
        temp_td.glm_ext_model = td_list{task_idx}.glm_ext_vel_model(in_x_dir,:,i_rep);
        temp_td.glm_handelbow_model = td_list{task_idx}.glm_handelbow_vel_model(in_x_dir,:,i_rep);
        eval_ext_x(:,i_rep) = evalModel(temp_td,struct('out_signals','LeftS1_FR','model_name','ext_model','model_type','glm','eval_metric','pr2'));
        eval_handelbow_x(:,i_rep) = evalModel(temp_td,struct('out_signals','LeftS1_FR','model_name','handelbow_model','model_type','glm','eval_metric','pr2'));
        
    end
    
    mean_ext_y = mean(eval_ext_y,2);
    mean_ext_x = mean(eval_ext_x,2);
    mean_handelbow_y = mean(eval_handelbow_y,2);
    mean_handelbow_x = mean(eval_handelbow_x,2);