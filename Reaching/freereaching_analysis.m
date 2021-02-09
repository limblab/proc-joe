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
        folderpath = 'D:\Lab\Data\DLC_videos\Han_20201204_rwFreeReach\neural-data\';
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
                fprintf('Removed %d percent of trials because of missing markers\n',sum(bad_points)/numel(bad_points)*100)
            end
            % remove trials where monkey's arm position is out of the
            % "workspace"
            td_list{i_td} = getExperimentPhase(td_list{i_td},task_list{i_td});

            
            
            % get robot height (z data of a hand marker). This is
            % meaningless during free reach...
            markername = 'hand3';
            dlc_idx = find((strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_z'])));
            robot_height(end+1) = mean(td_list{i_td}.dlc_pos(:,dlc_idx),'omitnan');
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
            'num_repeats',1,... % Raeed used 20 repeats.
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

    unit_idx = 25;
    td_idx = 1; % pick which trial data (task) to use

    window_plot = [-2,2]; % s
    num_trials_plot = 4;
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
    
    figure('Position',[145 558 1095 420]); hold on;
    % plot FR for unit
    x_data = (0:1:diff(window_idx(1,:)))*bin_size;
    encoderResults = encoderResults_cell{1};
    repeatnum = 1;
    % predict FR for unit in time window, then plot for each model
    for mdlnum = 1:numel(models_to_plot)
        td_list_smooth = getModel(td_list_smooth,encoderResults.glm_info{repeatnum,mdlnum});
    end
    for i_trial = 1:num_trials_plot
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
    subplot(1,4,1); 
    l=legend('Actual','Hand-only','Whole-arm');
    set(l,'box','off');
    xlabel('Time (s)');
    ylabel('FR (Hz)');
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
        for i_mdl = 1:length(models_to_plot) % 1 = 3D, 2 = 2D
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
                
            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
            axis square
            if monkeynum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel(spacenames(2))
            ylabel(spacenames(1))
        end
        legend(legend_data,models_to_plot)
    end
    suptitle('Within Pseudo-R^2 pairwise comparisons')
    
%% plot PDs across spaces (actual, as well as predicted by each model)
    sessionnum = 1;
    spacenames = {'RT3D','RT2D'};
    tasknames = {'Unknown','RW'};
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
