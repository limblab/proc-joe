%% load in TD_list
%% compute cross-correlation in sliding window for each trial data

    window_size = 10;
    window_step = 1;
    figure(); hold on;
    corr_list = {};
    for i_td = 1:numel(td_list)
        td_temp = td_list{i_td};
        dlc_idx = [find((strcmpi(td_temp.dlc_pos_names,['hand2','_x']))),...
                        find((strcmpi(td_temp.dlc_pos_names,['hand2','_y']))),...
                        find((strcmpi(td_temp.dlc_pos_names,['elbow1','_x']))),...
                        find((strcmpi(td_temp.dlc_pos_names,['elbow1','_y'])))];
                    
        window_starts = 1:window_step:(size(td_temp.dlc_pos,1)-window_size+1);
        
        corrs = zeros(size(window_starts,1),3);
        
        for i_win = 1:numel(window_starts)
            data = td_temp.dlc_vel(window_starts(i_win):window_starts(i_win)+window_size-1, dlc_idx);
            [r,lags] = xcorr(data(:,2),data(:,4),50,'biased');
            plot(lags,r);
            % get elbow-x hand-x corr
            corrs(i_win,1) = corr(data(:,3),data(:,1));
            % get elbow-y hand-y corr
            corrs(i_win,2) = corr(data(:,4),data(:,2));
            % get elbow-speed hand-speed corr
            corrs(i_win,3) = corr(sqrt(sum(data(:,1:2).^2,2)), sqrt(sum(data(:,3:4).^2,2)));
        end
        
        corr_list{i_td} = corrs;
        
    end


%% segment data into high correlation and low correlation segments. Make a TD mask
    % use hand,elbow y to start
td_high_corr = {};
td_low_corr = {};

    for i_td = 1:numel(td_list)
        corr_idx = 2;
        corr_use = corr_list{i_td}; % get 3D task
%         corr_cutoff = prctile(corr_use(:,corr_idx),50);
        corr_cutoff_high = prctile(corr_use(:,corr_idx),60);
        corr_cutoff_low = prctile(corr_use(:,corr_idx),40);
        
        corr_mask_high = corr_use(:,corr_idx) > corr_cutoff_high;
        corr_mask_high = [corr_mask_high; zeros(window_size-1,1)]; 

        corr_mask_low = corr_use(:,corr_idx) < corr_cutoff_low;
        corr_mask_low = [corr_mask_low; zeros(window_size-1,1)]; 
        
        td_high_corr_temp = td_list{i_td};
        td_low_corr_temp = td_list{i_td};

        % remove data points based on correlation and make data lengths the
        % same across files
        
        field_len = length(td_high_corr_temp.dlc_pos);
        [desired_field_len,desired_idx] = min([sum(corr_mask_high), sum(corr_mask_low)]);
        td_fieldnames = fieldnames(td_high_corr_temp);

        keep_idx = datasample(1:1:max([sum(corr_mask_high), sum(corr_mask_low)]),desired_field_len,'Replace',false);

        
        for i_field = 1:numel(td_fieldnames)
            if(length(td_high_corr_temp.(td_fieldnames{i_field})) == field_len)
                td_high_corr_temp.(td_fieldnames{i_field}) = td_high_corr_temp.(td_fieldnames{i_field})(corr_mask_high==1,:);
                if(desired_idx == 2)
                    td_high_corr_temp.(td_fieldnames{i_field}) = td_high_corr_temp.(td_fieldnames{i_field})(keep_idx,:);
                end
            end
        end    

        field_len = length(td_low_corr_temp.dlc_pos);
        td_fieldnames = fieldnames(td_low_corr_temp);

        for i_field = 1:numel(td_fieldnames)
            if(length(td_low_corr_temp.(td_fieldnames{i_field})) == field_len)
                td_low_corr_temp.(td_fieldnames{i_field}) = td_low_corr_temp.(td_fieldnames{i_field})(corr_mask_low==1,:);
                if(desired_idx == 1)
                    td_low_corr_temp.(td_fieldnames{i_field}) = td_low_corr_temp.(td_fieldnames{i_field})(keep_idx,:);
                end
            end
        end 
               
        
        td_high_corr{end+1} = td_high_corr_temp;
        td_low_corr{end+1} = td_low_corr_temp;
    end
    
    
    
    
%% Get encoding models for high and low correlation segments

    arrayname = 'LeftS1';
    monkey_names = {'Han'};
    included_models = {'ext','handelbow'}; % models to calculate encoders for
    models_to_plot = included_models; % main models of the paper
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


    encoderResults_cell = {};
    for i=1:2
        switch i
            case 1
                td_use = td_high_corr;
            case 2 
                td_use = td_low_corr;
        end

        encoderResults_cell{i} = reachingEncoders(td_use,task_list,0,struct(...
            'model_aliases',{included_models},...
            'arrayname',arrayname,...
            'num_tuning_bins',16,...
            'crossval_lookup',[],...
            'get_tuning_curves',true,...
            'num_repeats',5,... % Raeed used 20 repeats.
            'num_folds',5));

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
    
%% make plots for each file (high and low correlation in this case)
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
    spacenames = {'RT3D','RT2D'};

    legend_data = [];

%     f.Name = [encoderResults.crossEval.monkey{1} ,'_', encoderResults.crossEval.date{1},'_encoderModelComparison'];
    for monkeynum = 1:length(monkey_names)
        for spacenum = 1:2 % 1 = 3D, 2 = 2D
            f=figure();
            % set subplot
            subplot(1,length(monkey_names),monkeynum)
            plot([-1 1],[-1 1],'k--','linewidth',0.5)
            hold on
            plot([0 0],[-1 1],'k-','linewidth',0.5)
            plot([-1 1],[0 0],'k-','linewidth',0.5)

            avg_pR2_high = neuronAverage(model_eval{monkeynum,1},struct('keycols','signalID','do_ci',false));
            avg_pR2_low = neuronAverage(model_eval{monkeynum,2},struct('keycols','signalID','do_ci',false));   
            model_extension = ['_space',num2str(spacenum),'_within_eval'];

            legend_data(end+1) = scatter(...
                avg_pR2_low.(strcat(model_pairs{1,2},model_extension)),... % switch between models in model_pairs
                avg_pR2_high.(strcat(model_pairs{1,2},model_extension)),...
                [],getColorFromList(1,spacenum+1),'filled');


            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
            axis square
            if monkeynum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
%             xlabel(sprintf('%s pR2',getModelTitles(model_pairs{1,1})))
%             ylabel(sprintf('%s pR2',getModelTitles(model_pairs{1,2})))
            xlabel('Low pR2');
            ylabel('High pR2');
        end
    end
    formatForLee(gcf);
    set(gca,'fontsize',14);

