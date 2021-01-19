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
        folderpath = 'D:\Lab\Data\DLC_videos\Han_20201223_rwTwoPlanes\neural-data';
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
    robot_height_all = {};
    for filenum = 1:length(filenames)
        % Load data
        load(fullfile(file_info(filenum).folder,file_info(filenum).name));
        keep_mask = ones(size(td_list));
        robot_height = [];
        for i_td = 1:numel(td_list) % each entry in td_list is a different trial_data for a different experiment (free reach vs 2D random walk for example)            
            % resample trial data to appropriate bin size
            if(strcmpi(task_list{i_td},'RT3D')==1)
                keep_mask(i_td) = 0;
            else
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
                    % get marker velocity
                    td_list{i_td} = getDifferential(td_list{i_td},struct('signals','dlc_pos','alias','dlc_vel'));
                    % remove time points where dlc tracking is bad
                    bad_points = any(isnan(td_list{i_td}.dlc_pos),2) | any(isnan(td_list{i_td}.dlc_vel),2);
                    td_names = fieldnames(td_list{i_td});
                    for i_name = 1:numel(td_names)
                        if(size(bad_points,1) == size(td_list{i_td}.(td_names{i_name}),1))
                            td_list{i_td}.(td_names{i_name})(bad_points,:) = [];
                        end
                    end
                    fprintf('Removed %d percent of trials because of missing markers\n',sum(bad_points)/numel(bad_points))
                end
                % remove trials where monkey's arm position is out of the
                % "workspace"
                td_list{i_td} = getExperimentPhase(td_list{i_td},task_list{i_td});
                
                % get robot height (z data of a hand marker)
                markername = 'hand3';
                dlc_idx = find((strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_z'])));
                robot_height(end+1) = mean(td_list{i_td}.dlc_pos(:,dlc_idx));
            end
        end
        td_all{filenum} = td_list(keep_mask==1);
        robot_height_all{filenum} = robot_height;
        clear td_list;
    end

    
    
%% plot handle vs hand position
    markername = 'hand3';
    for filenum = 1:numel(td_all)
        td_list = td_all{filenum};
        for i_td = 1:numel(td_list)
            if(strcmpi(task_list{i_td},'RT'))
                figure();
                dlc_idx = [find((strcmpi(td_list{1}.dlc_pos_names,[markername,'_x']))),find((strcmpi(td_list{1}.dlc_pos_names,[markername,'_y'])))];
                plot(td_list{i_td}.pos(:,2),td_list{i_td}.dlc_pos(:,dlc_idx(1)),'k');
                hold on;
                
            end
        end
    end

    
%% Loop through files to cross-validate encoders
    fprintf('Starting analysis of %d files. Warning: cross-validation will take a long time\n',length(td_all))
    encoderResults_cell = cell(length(td_all),1);
    for filenum = 1:length(td_all)
        % Load data
        td_list = td_all{filenum};
        
        % bin data at 50ms
        for i_td = 1:numel(td_list)
            td_list{i_td} = binTD(td_list{i_td},0.05/td_list{i_td}(1).bin_size);
        end
        
        % Get encoding models
        encoderResults_cell{filenum} = twoPlanesEncoders(td_list,robot_height_all{filenum},struct(...
            'model_aliases',{included_models},...
            'arrayname',arrayname,...
            'num_tuning_bins',16,...
            'crossval_lookup',[],...
            'get_tuning_curves',true,...
            'num_repeats',2,...
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
        [space_eval_cell,space_eval_within_cell] = deal(cell(2,length(included_models)));
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
            end
        end
        model_eval{monkey_idx,session_ctr(monkey_idx)} = horzcat(...
            model_eval{monkey_idx,session_ctr(monkey_idx)},...
            model_eval_cell{:},...
            space_eval_cell{:},...
            space_eval_within_cell{:});

        % We already have tuning table in crossTuning... just extract the models we want
        model_tuning{monkey_idx,session_ctr(monkey_idx)} = encoderResults.crossTuning(:,...
            contains(encoderResults.crossTuning.Properties.VariableDescriptions,'meta') |...
            strcmpi(encoderResults.crossTuning.Properties.VariableNames,'bins'));
        model_tuning_cell = cell(1,length(included_models)+1);
        for modelnum = 1:length(included_models)
            model_tuning_cell{modelnum} = table(...
                encoderResults.crossTuning.(sprintf('glm_%s_model_velCurve',included_models{modelnum})),...
                encoderResults.crossTuning.(sprintf('glm_%s_model_velPD',included_models{modelnum})),...
                'VariableNames',strcat(included_models(modelnum),{'_velCurve','_velPD'}));
            model_tuning_cell{modelnum}.Properties.VariableDescriptions = {'linear','circular'};
        end
        model_tuning_cell{end} = table(...
            encoderResults.crossTuning.([arrayname,'_FR_velCurve']),...
            encoderResults.crossTuning.([arrayname,'_FR_velPD']),...
            'VariableNames',strcat([arrayname,'_FR'],{'_velCurve','_velPD'}));
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
                encoderResults,struct('model_aliases',{included_models},'neural_signal',[arrayname,'_FR']));
            tuned_neurons{monkey_idx,session_ctr(monkey_idx)} = encoderResults.tunedNeurons;
        else
            warning('No tuned neurons field found!')
        end

        % output a counter
        fprintf('Processed file %d of %d at time %f\n',filenum,length(encoderResults_cell),toc(fileclock))
    end
    
 
%% Get example tuning curves for all models
    for monkeynum = 1:length(monkey_names)
        for sessionnum = 1:session_ctr(monkeynum)
            % Plot out tuning curves for tuned neurons
            f = figure('defaultaxesfontsize',18);
            % plot a max of 7 neurons
            num_neurons = min(7,size(tuned_neurons{monkeynum,sessionnum},1));
            neurons_to_plot = tuned_neurons{monkeynum,sessionnum}(...
                randperm(size(tuned_neurons{monkeynum,sessionnum},1),num_neurons),...
                :);
            for neuronnum = 1:num_neurons
                % figure out maxFR over both workspaces
                [~,temp_table] = getNTidx(model_tuning{monkeynum,sessionnum},...
                    'signalID',neurons_to_plot(neuronnum,:));
                minFR = floor(min(min(temp_table{:,strcat([models_to_plot {[arrayname,'_FR']}],'_velCurve')})));
                maxFR = ceil(max(max(temp_table{:,strcat([models_to_plot {[arrayname,'_FR']}],'_velCurve')})));

                % go plot each workspace
                for spacenum = 1:size(cond_colors,1)
                    % get tuning table specifically for this neuron and workspace
                    [~,temp_table] = getNTidx(model_tuning{monkeynum,sessionnum},...
                        'signalID',neurons_to_plot(neuronnum,:),...
                        'spaceNum',spacenum);

                    % plot out actual tuning curves
                    subplot(length(models_to_plot)+1,num_neurons,neuronnum)
                    plotTuning(temp_table,...
                        struct('maxFR',maxFR,...
                            'minFR',minFR,...
                            'unroll',true,...
                            'color',cond_colors(spacenum,:),...
                            'curve_colname',sprintf('%s_velCurve',[arrayname,'_FR']),...
                            'pd_colname',sprintf('%s_velPD',[arrayname,'_FR']),...
                            'plot_ci',false))
                    title(strcat('Neuron ',num2str(neurons_to_plot(neuronnum,:))))

                    % plot out modeled tuning curves
                    for modelnum = 1:length(models_to_plot)
                        subplot(length(models_to_plot)+1,num_neurons,modelnum*num_neurons+neuronnum)
                        plotTuning(temp_table,...
                            struct('maxFR',maxFR,...
                                'minFR',minFR,...
                                'unroll',true,...
                                'color',cond_colors(spacenum,:),...
                                'curve_colname',sprintf('%s_velCurve',models_to_plot{modelnum}),...
                                'pd_colname',sprintf('%s_velPD',models_to_plot{modelnum}),...
                                'plot_ci',false))
                        title(getModelTitles(models_to_plot{modelnum}),'interpreter','none')
                    end
                end
            end
            set(gcf,'renderer','Painters')
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

    % make the pairwise comparison scatter plot
    figure
    for monkeynum = 1:length(monkey_names)
        for pairnum = 1:size(model_pairs,1)
            % set subplot
            subplot(size(model_pairs,1),length(monkey_names),...
                (pairnum-1)*length(monkey_names)+monkeynum)
            plot([-1 1],[-1 1],'k--','linewidth',0.5)
            hold on
            plot([0 0],[-1 1],'k-','linewidth',0.5)
            plot([-1 1],[0 0],'k-','linewidth',0.5)
            for sessionnum = 1:session_ctr(monkeynum)
                avg_pR2 = neuronAverage(model_eval{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
                % scatter filled circles if there's a winner, empty circles if not
                no_winner =  cellfun(@isempty,pr2_winners{monkeynum,sessionnum}(pairnum,:));
                scatter(...
                    avg_pR2.(strcat(model_pairs{pairnum,1},'_eval'))(no_winner),...
                    avg_pR2.(strcat(model_pairs{pairnum,2},'_eval'))(no_winner),...
                    [],session_colors(sessionnum,:))
                scatter(...
                    avg_pR2.(strcat(model_pairs{pairnum,1},'_eval'))(~no_winner),...
                    avg_pR2.(strcat(model_pairs{pairnum,2},'_eval'))(~no_winner),...
                    [],session_colors(sessionnum,:),'filled')
            end
            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
            axis square
            if monkeynum ~= 1 || pairnum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,1})))
            ylabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,2})))
        end
    end
    suptitle('Pseudo-R^2 pairwise comparisons')

    % show scatter plot for hand/elbow pR2 within condition vs against condition
    for modelnum = 1:length(models_to_plot)
        figure
        for monkeynum = 1:length(monkey_names)
            for spacenum = 1:2
                % set subplot
                subplot(2,length(monkey_names),(spacenum-1)*length(monkey_names)+monkeynum)

                % plot lines
                plot([-1 1],[-1 1],'k--','linewidth',0.5)
                hold on
                plot([0 0],[-1 1],'k-','linewidth',0.5)
                plot([-1 1],[0 0],'k-','linewidth',0.5)

                % plot out each session
                for sessionnum = 1:session_ctr(monkeynum)
                    avg_pR2 = neuronAverage(model_eval{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
                    scatter(...
                        avg_pR2.(sprintf('%s_space%d_eval',models_to_plot{modelnum},spacenum)),...
                        avg_pR2.(sprintf('%s_space%d_within_eval',models_to_plot{modelnum},spacenum)),...
                        [],session_colors(sessionnum,:),'filled')
                end
                % make axes pretty
                set(gca,'box','off','tickdir','out',...
                    'xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
                axis square
                if monkeynum ~= 1 || spacenum ~= 1
                    set(gca,'box','off','tickdir','out',...
                        'xtick',[],'ytick',[])
                end
                xlabel(sprintf('%s trained across pR2',getModelTitles(models_to_plot{modelnum})))
                ylabel(sprintf('%s trained within pR2',getModelTitles(models_to_plot{modelnum})))
                title(sprintf('Workspace %d',spacenum))
            end
        end
        suptitle('Full pR^2 vs within condition pR^2')
    end

%% Tuning curve shape comparison
    % find winners of tuning corr
    tuning_corr_winners = cell(length(monkey_names),size(session_colors,1));
    for monkeynum = 1:length(monkey_names)
        for sessionnum = 1:session_ctr(monkeynum)
            [tuning_corr_winners{monkeynum,sessionnum},model_pairs] = compareEncoderMetrics(...
                    tuning_corr{monkeynum,sessionnum},struct(...
                        'bonferroni_correction',6,...
                        'models',{models_to_plot},...
                        'model_pairs',{{'ext','handelbow'}},...
                        'postfix','_tuningCorr'));
        end
    end

    % figure out how many neurons the hand-based models could beat either of the whole-arm models
    all_tuning_corr_winners = horzcat(tuning_corr_winners{:});
    ext_winners = sum(strcmpi(all_tuning_corr_winners,'ext'),2);
    handelbow_winners = sum(strcmpi(all_tuning_corr_winners,'handelbow'),2);
    fprintf('tuning correlation winners -- hand-only: %d, whole-arm: %d\n',ext_winners,handelbow_winners)

    % make the pairwise comparison scatter plot
    figure
    for monkeynum = 1:length(monkey_names)
        for pairnum = 1:size(model_pairs,1)
            % set subplot
            subplot(size(model_pairs,1),length(monkey_names),...
                (pairnum-1)*length(monkey_names)+monkeynum)
            plot([-1 1],[-1 1],'k--','linewidth',0.5)
            hold on
            plot([0 0],[-1 1],'k-','linewidth',0.5)
            plot([-1 1],[0 0],'k-','linewidth',0.5)
            for sessionnum = 1:session_ctr(monkeynum)
                avg_corr = neuronAverage(tuning_corr{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
                % scatter filled circles if there's a winner, empty circles if not
                no_winner =  cellfun(@isempty,tuning_corr_winners{monkeynum,sessionnum}(pairnum,:));
                scatter(...
                    avg_corr.(strcat(model_pairs{pairnum,1},'_tuningCorr'))(no_winner),...
                    avg_corr.(strcat(model_pairs{pairnum,2},'_tuningCorr'))(no_winner),...
                    [],session_colors(sessionnum,:))
                scatter(...
                    avg_corr.(strcat(model_pairs{pairnum,1},'_tuningCorr'))(~no_winner),...
                    avg_corr.(strcat(model_pairs{pairnum,2},'_tuningCorr'))(~no_winner),...
                    [],session_colors(sessionnum,:),'filled')
            end
            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-0.1 1],'ylim',[-0.1 1],...
                'xtick',0:0.5:1,'ytick',0:0.5:1)
            axis square
            if monkeynum ~= 1 || pairnum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel(sprintf('%s',getModelTitles(model_pairs{pairnum,1})))
            ylabel(sprintf('%s',getModelTitles(model_pairs{pairnum,2})))
        end
    end
    suptitle('Tuning correlation pairwise comparisons')

%% PD shifts over all monkeys
    file_shifts = cell(length(encoderResults_cell),length(models_to_plot)); % shift tables for each model in each file
    for filenum = 1:length(encoderResults_cell)
        % load data
        encoderResults = encoderResults_cell{filenum};

        shift_tables = calculatePDShiftTables(encoderResults,[strcat('glm_',models_to_plot,'_model'), [arrayname,'_FR']]);
        mean_shifts = cell(length(models_to_plot),1);
        for modelnum = 1:length(models_to_plot)+1
            mean_shifts{modelnum} = neuronAverage(shift_tables{modelnum},struct(...
                'keycols',{{'monkey','date','task','signalID'}}));
            [~,file_shifts{filenum,modelnum}] = getNTidx(mean_shifts{modelnum},'signalID',encoderResults.tunedNeurons);
        end
    end

    allFileShifts_real = vertcat(file_shifts{:,end});

    % Make histograms and scatters
    total_hists = figure('defaultaxesfontsize',18);
    scatters = figure('defaultaxesfontsize',18);
    for monkeynum = 1:length(monkey_names)
        % get monkey specific session dates
        [~,monkey_shifts_real] = getNTidx(allFileShifts_real,'monkey',monkey_names{monkeynum});
        session_dates = unique(monkey_shifts_real.date);

        % make histogram combining sessions
            % first figure out ylim
            % ylim_high = 10*floor(height(monkey_shifts_real)/20);
            ylim_high = 40;
            % actual PD shift histogram
            figure(total_hists)
            subplot(length(models_to_plot)+1,length(monkey_names),monkeynum)
            h = histogram(gca,monkey_shifts_real.velPD*180/pi,'BinWidth',10,'DisplayStyle','stair');
            set(h,'facecolor','none','edgecolor',ones(1,3)*0.5)
            set(gca,...
                'box','off','tickdir','out',...
                'xlim',[-180 180],'xtick',[-180 0 180],...
                'ylim',[0 ylim_high],'ytick',[0 ylim_high/2 ylim_high])
            title(sprintf('Monkey %s',monkey_names{monkeynum}))
            if monkeynum == 1
                ylabel('Actual PD Shift')
            end
            for modelnum = 1:length(models_to_plot)
                allFileShifts_model = vertcat(file_shifts{:,modelnum});
                [~,monkey_shifts_model] = getNTidx(allFileShifts_model,'monkey',monkey_names{monkeynum});

                % modeled PD shift histogram
                subplot(length(models_to_plot)+1,length(monkey_names),modelnum*length(monkey_names)+monkeynum)
                h = histogram(gca,monkey_shifts_model.velPD*180/pi,'BinWidth',10,'DisplayStyle','stair');
                set(h,'facecolor','none','edgecolor',ones(1,3)*0.5)
                set(gca,...
                    'box','off','tickdir','out',...
                    'xlim',[-180 180],'xtick',[-180 0 180],...
                    'ylim',[0 ylim_high],'ytick',[0 ylim_high/2 ylim_high])
                if monkeynum == 1
                    ylabel(sprintf('%s modeled PD shift',getModelTitles(models_to_plot{modelnum})))
                end
            end

        % make scatter plots separating sessions
            for sessionnum = 1:length(session_dates)
                % get real shifts for this session
                [~,session_shifts_real] = getNTidx(allFileShifts_real,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});

                % now the models
                for modelnum = 1:length(models_to_plot)
                    % get the modeled shifts for this session
                    allFileShifts_model = vertcat(file_shifts{:,modelnum});
                    [~,session_shifts_model] = getNTidx(allFileShifts_model,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});

                    % scatter plots
                    figure(scatters)
                    subplot(length(models_to_plot),length(monkey_names),(modelnum-1)*length(monkey_names)+monkeynum)
                    scatter(...
                        180/pi*session_shifts_real.velPD,...
                        180/pi*session_shifts_model.velPD,...
                        50,session_colors(sessionnum,:),'filled')
                    hold on
                end
            end

            % make plots pretty
            for modelnum = 1:length(models_to_plot)
                % scatter plots
                figure(scatters)
                subplot(length(models_to_plot),length(monkey_names),(modelnum-1)*length(monkey_names)+monkeynum)
                plot([-180 180],[0 0],'-k','linewidth',2)
                plot([0 0],[-180 180],'-k','linewidth',2)
                plot([-180 180],[-180 180],'--k','linewidth',2)
                axis equal
                set(gca,...
                    'box','off','tickdir','out',...
                    'xtick',[-180 180],'ytick',[-180 180],...
                    'xlim',[-180 180],'ylim',[-180 180])
                % labels
                if modelnum == length(models_to_plot)
                    xlabel('Actual PD Shift')
                end
                if modelnum == 1
                    title(sprintf('Monkey %s',monkey_names{monkeynum}))
                end
                if monkeynum == 1
                    ylabel(...
                        {sprintf('%s model',getModelTitles(models_to_plot{modelnum}));'Modeled PD Shift'},...
                        'interpreter','none')
                end
            end
    end
    figure(total_hists)
    suptitle('PD shift histograms')
    figure(scatters)
    suptitle('PD shift scatter plots')

    % PD shift VAF dotplots
        % find winners of PD shift
        shift_vaf_winners = cell(length(monkey_names),size(session_colors,1));
        for monkeynum = 1:length(monkey_names)
            for sessionnum = 1:session_ctr(monkeynum)
                [shift_vaf_winners{monkeynum,sessionnum},model_pairs] = compareEncoderMetrics(...
                        shift_vaf{monkeynum,sessionnum},struct(...
                            'bonferroni_correction',6,...
                            'models',{models_to_plot},...
                            'model_pairs',{{'ext','handelbow'}},...
                            'postfix','_vaf'));
            end
        end

        % Find session winners of PD shift
        shift_vaf_session_winners = cell(length(monkey_names),size(session_colors,1));
        for monkeynum = 1:length(monkey_names)
            for sessionnum = 1:session_ctr(monkeynum)
                shift_vaf_session = neuronAverage(shift_vaf{monkeynum,sessionnum},struct(...
                    'keycols',{{'monkey','date','task','crossvalID'}},...
                    'do_ci',false));
                temp_tab = table(ones(height(shift_vaf_session),1),'VariableNames',{'signalID'});
                shift_vaf_session = horzcat(temp_tab,shift_vaf_session);
                [shift_vaf_session_winners{monkeynum,sessionnum},model_pairs] = compareEncoderMetrics(...
                        shift_vaf_session,struct(...
                            'bonferroni_correction',6,...
                            'models',{models_to_plot},...
                            'model_pairs',{{'ext','handelbow'}},...
                            'postfix','_vaf'));
            end
        end

        % find winners
        all_shift_vaf_winners = horzcat(shift_vaf_winners{:});
        ext_winners = sum(strcmpi(all_shift_vaf_winners,'ext'),2);
        handelbow_winners = sum(strcmpi(all_shift_vaf_winners,'handelbow'),2);
        fprintf('shift VAF winners -- hand-only: %d, whole-arm: %d\n',ext_winners,handelbow_winners)

        % session winners
        all_shift_vaf_session_winners = horzcat(shift_vaf_session_winners{:});
        ext_winners = sum(strcmpi(all_shift_vaf_session_winners,'ext'),2);
        handelbow_winners = sum(strcmpi(all_shift_vaf_session_winners,'handelbow'),2);
        fprintf('shift VAF session winners -- hand-only: %d, whole-arm: %d\n',ext_winners,handelbow_winners)

        % get average shift vaf over all neurons
        shift_vaf_all = vertcat(shift_vaf{:});
        avg_shift_vaf_all = neuronAverage(shift_vaf_all,struct(...
            'keycols',{{'task'}},...
            'do_ci',false));

        % plot session averages with CI bars
        figure('defaultaxesfontsize',18)
        alpha = 0.05;
        model_x = (2:3:((length(models_to_plot)-1)*3+2))/10;
        for monkeynum = 1:length(monkey_names)
            subplot(length(monkey_names),1,monkeynum)
            for sessionnum = 1:session_ctr(monkeynum)
                % plot session average
                avg_shift_vaf = neuronAverage(shift_vaf{monkeynum,sessionnum},...
                    struct('keycols',{{'monkey','date','task','crossvalID'}},'do_ci',false));

                % estimate error bars
                [~,cols] = ismember(strcat(models_to_plot,'_vaf'),avg_shift_vaf.Properties.VariableNames);
                num_repeats = double(max(shift_vaf{monkeynum,sessionnum}.crossvalID(:,1)));
                num_folds = double(max(shift_vaf{monkeynum,sessionnum}.crossvalID(:,2)));
                crossval_correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
                yvals = mean(avg_shift_vaf{:,cols});
                var_vaf = var(avg_shift_vaf{:,cols});
                upp = tinv(1-alpha/2,num_folds*num_repeats-1);
                low = tinv(alpha/2,num_folds*num_repeats-1);
                CI_lo = yvals + low * sqrt(crossval_correction*var_vaf);
                CI_hi = yvals + upp * sqrt(crossval_correction*var_vaf);
                
                % plot dots and lines
                session_jitter = 0.01*(sessionnum-(session_ctr(monkeynum)+1)/2);
                plot(repmat(model_x+session_jitter,2,1),[CI_lo;CI_hi],'-','linewidth',2,'color',session_colors(sessionnum,:))
                hold on
                scatter(model_x(:)+session_jitter,yvals(:),50,session_colors(sessionnum,:),'filled')
            end
            ylabel('PD Shift Circular VAF')
            set(gca,'box','off','tickdir','out',...
                'xlim',[model_x(1)-0.2 model_x(end)+0.2],'xtick',model_x,'xticklabel',getModelTitles(models_to_plot),...
                'ylim',[0 1],'ytick',[0 1])
        end

