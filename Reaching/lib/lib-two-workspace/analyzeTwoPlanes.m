%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [crossEval, crossTuning, crossvalLookup] = analyzeTwoPlanes(trial_data,params)
% 
% For multiworkspace files, with dl and pm workspaces:
%   * Fits three different coordinate frame models to data from both workspaces
%   * Predicts firing rates from models on test data in both workspaces
%   * Compare predicted firing rates to actual firing rates in test data
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .arrayname  : which signals to calculate PDs for
%                           default: 'S1'
%       .model_type     :   type of model to fit
%                           default; 'glm'
%       .glm_distribution : distribution to use for GLM
%                           default: 'poisson'
%       .model_eval_metric : Evaluation metric for models
%                           default: 'pr2'
%       .num_musc_pcs   :   Number of muscle principle components to use
%                           default: 5
%       .num_repeats    : # crossval repeats to use
%                           default: 20
%       .num_folds    : # crossval folds to use
%                           default: 5
%       .crossval_lookup    : lookup table to pick specific trials for each crossval
%                           (columns for crossvalID and trialID)
%                           default: []
%       .verbose      :     Print diagnostic information
%                           default: true
%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [crossEval, crossTuning, crossvalLookup] = analyzeTwoPlanes(td_list,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up
    % default parameters
    num_folds = 5;
    num_repeats = 20;
    crossval_lookup = [];
    verbose = true;
    if nargin > 1, assignParams(who,params); end % overwrite parameters

    glm_params = params.glm_params;
%% Compile training and test sets
    % inialize temporary eval holders
    [repeatEval,repeatTuning,repeatCrossvalLookup] = deal(cell(num_repeats,1));

    % extract td_high and td_low
    td_high = td_list{1};
    td_low = td_list{2};

    % loop over num repeats
    if verbose
        repeat_timer = tic;
        fprintf('Starting %dx%d-fold crossvalidation at time %f\n',num_repeats,num_folds,toc(repeat_timer));
    end
    for repeatctr = 1:num_repeats
        % get fold indices
        indices = crossvalind('Kfold',length(td_high.pos),num_folds);

        % initialize temporary fold evaluation structure
        [foldEval,foldTuning,foldCrossvalLookup] = deal(cell(num_folds,1));

        % loop over number of folds
        if verbose
            fold_timer = tic;
        end
        for foldctr = 1:num_folds
            % Get test and training indices for this fold
            % check if there's a crossval lookup table
            td_test_low = []; td_test_high = [];
            if isempty(crossval_lookup)
                test_idx = (indices==foldctr);
                train_idx = ~test_idx;
                % split into td_train and td_test, only use in_signals and
                % out_signals. Also grab handle velocity for td_test so that we
                % can compute PDs using handle kinematics
                for modelnum = 1:numel(glm_params)
                    for i_name = 1:size(glm_params{modelnum}.in_signals,2)
                        for i_idx = 1:numel(glm_params{modelnum}.in_signals{i_name,2})
                            td_train.(glm_params{modelnum}.in_signals{i_name,1}) = ...
                                [td_high.(glm_params{modelnum}.in_signals{i_name,1})(train_idx,:); td_low.(glm_params{modelnum}.in_signals{i_name,1})(train_idx,:)];
                            td_test_high.(glm_params{modelnum}.in_signals{i_name,1}) = td_high.(glm_params{modelnum}.in_signals{i_name,1})(test_idx,:);
                            td_test_low.(glm_params{modelnum}.in_signals{i_name,1}) = td_low.(glm_params{modelnum}.in_signals{i_name,1})(test_idx,:);
                            % meta info....
                            td_train.monkey = td_high.monkey;
                            td_test_high.monkey = td_high.monkey;
                            td_test_low.monkey = td_low.monkey;
            
                            td_train.task = td_high.task;
                            td_test_high.task = td_high.task;
                            td_test_low.task = td_low.task;

                            td_train.date = td_high.date;
                            td_test_high.date = td_high.date;
                            td_test_low.date = td_low.date;
                        end
                    end
                end
                td_test_high.vel = td_high.vel(test_idx,:);
                td_test_low.vel = td_low.vel(test_idx,:);
                
                % get out signals
                td_train.(glm_params{1}.out_signals) = ...
                    [td_high.(glm_params{1}.out_signals)(train_idx,:); td_low.(glm_params{1}.out_signals)(train_idx,:)];
                td_test_high.(glm_params{1}.out_signals) = td_high.(glm_params{1}.out_signals)(test_idx,:); 
                td_test_low.(glm_params{1}.out_signals) = td_low.(glm_params{1}.out_signals)(test_idx,:);

                % set bin_size as well
                td_train.bin_size = td_high.bin_size;
                td_test_high.bin_size = td_low.bin_size;
                td_test_low.bin_size = td_low.bin_size;
                
                % set spacenum (1=high, 2=low)
                train_spacenum = [ones(sum(train_idx),1);1+ones(sum(train_idx),1)];
            else
                error('crossval lookup is not empty. This is not implemented');
            end

            % analyze fold to get model evaluations
            if exist('params','var')
                params.crossvalID = uint16([repeatctr foldctr]);
            else
                params = struct('crossvalID',uint16([repeatctr foldctr]));
            end
            [foldEval{foldctr},foldTuning{foldctr}] = analyzeFold(td_train,{td_test_high,td_test_low},train_spacenum,params);

            % get test trialIDs
            trialID = cat(1,find(indices==foldctr)',find(indices==foldctr)');
            spaceNum = zeros(size(trialID)) + [1;2]; % 1=high,2=low
            foldCrossvalLookup{foldctr} = table(...
                repmat(params.crossvalID,length(trialID),1),trialID',spaceNum',...
                'VariableNames',{'crossvalID','trialID','spaceNum'});

            if verbose
                fprintf('\tEvaluated fold %d of %d at time %f\n',foldctr,num_folds,toc(fold_timer));
            end
        end

        % put fold outputs into larger table
        repeatEval{repeatctr} = vertcat(foldEval{:});
        repeatTuning{repeatctr} = vertcat(foldTuning{:});
        repeatCrossvalLookup{repeatctr} = vertcat(foldCrossvalLookup{:});

        if verbose
            fprintf('Evaluated repeat %d of %d at time %f\n',repeatctr,num_repeats,toc(repeat_timer));
        end
    end

    % put all evals together
    crossEval = vertcat(repeatEval{:});
    crossTuning = vertcat(repeatTuning{:});
    crossvalLookup = vertcat(repeatCrossvalLookup{:});

%% Diagnostics...
    % [foldEval,foldTuning] = analyzeFold(td_train,td_test);
    % musc_err = minusPi2Pi(foldEval.glm_musc_model_velPDShift-foldEval.S1_FR_velPDShift);
    % ext_err = minusPi2Pi(foldEval.glm_ext_model_velPDShift-foldEval.S1_FR_velPDShift);
    % ego_err = minusPi2Pi(foldEval.glm_ego_model_velPDShift-foldEval.S1_FR_velPDShift);

    % err_frac_musc = circ_var(musc_err)/circ_var(foldEval.S1_FR_velPDShift);
    % err_frac_ext = circ_var(ext_err)/circ_var(foldEval.S1_FR_velPDShift);
    % err_frac_ego = circ_var(ego_err)/circ_var(foldEval.S1_FR_velPDShift);

    % figure
    % scatter(foldEval.S1_FR_velPDShift,foldEval.glm_musc_model_velPDShift,[],'b')
    % hold on
    % scatter(foldEval.S1_FR_velPDShift,foldEval.glm_ext_model_velPDShift,[],'r')
    % scatter(foldEval.S1_FR_velPDShift,foldEval.glm_ego_model_velPDShift,[],'g')
    % plot([-pi pi],[-pi pi],'--k','linewidth',2)
    % axis equal

    % figure
    % polar(musc_err,ones(size(musc_err)),'bo')
    % hold on
    % polar([0 circ_mean(musc_err)],[0 circ_r(musc_err)],'b-')
    % set(gca,'xlim',[-1 1],'ylim',[-1 1])
    % axis equal
    % figure
    % polar(ext_err,ones(size(ext_err)),'ro')
    % hold on
    % polar([0 circ_mean(ext_err)],[0 circ_r(ext_err)],'r-')
    % set(gca,'xlim',[-1 1],'ylim',[-1 1])
    % axis equal
    % figure
    % polar(ego_err,ones(size(ego_err)),'go')
    % hold on
    % polar([0 circ_mean(ego_err)],[0 circ_r(ego_err)],'g-')
    % set(gca,'xlim',[-1 1],'ylim',[-1 1])
    % axis equal

    % figure
    % subplot(2,1,1)
    % scatter(cos(foldEval.S1_FR_velPDShift),cos(foldEval.glm_musc_model_velPDShift),[],'b')
    % hold on
    % plot([-1 1],[-1 1],'--k','linewidth',2)
    % axis equal
    % subplot(2,1,2)
    % scatter(sin(foldEval.S1_FR_velPDShift),sin(foldEval.glm_musc_model_velPDShift),[],'b')
    % hold on
    % plot([-1 1],[-1 1],'--k','linewidth',2)
    % axis equal

    % figure
    % subplot(2,1,1)
    % scatter(cos(foldEval.S1_FR_velPDShift),cos(foldEval.glm_ext_model_velPDShift),[],'r')
    % hold on
    % plot([-1 1],[-1 1],'--k','linewidth',2)
    % axis equal
    % subplot(2,1,2)
    % scatter(sin(foldEval.S1_FR_velPDShift),sin(foldEval.glm_ext_model_velPDShift),[],'r')
    % hold on
    % plot([-1 1],[-1 1],'--k','linewidth',2)
    % axis equal

    % figure
    % subplot(2,1,1)
    % scatter(cos(foldEval.S1_FR_velPDShift),cos(foldEval.glm_ego_model_velPDShift),[],'g')
    % hold on
    % plot([-1 1],[-1 1],'--k','linewidth',2)
    % axis equal
    % subplot(2,1,2)
    % scatter(sin(foldEval.S1_FR_velPDShift),sin(foldEval.glm_ego_model_velPDShift),[],'g')
    % hold on
    % plot([-1 1],[-1 1],'--k','linewidth',2)
    % axis equal

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [foldEval,foldTuning] = analyzeFold(td_train,td_test,train_spacenum,params)
% ANALYZEFOLD analyze single fold of cross-validation set and return NeuronTable structure
%   with evaluation information for this fold
% 
% Inputs:
%   td_train - trial data structure with trials to be used for training
%   td_test - cell array of trial_data structures to be used for testing two workspaces
%   params - parameters struct
%       .model_eval_metric - metric for model evaluation (default: pR2)
%       .glm_params - cell array of paramter structs to fit glm models with getModel
%       .model_names - names of models to check classical tuning for
%           (including actual neural signal)
%       .num_tuning_bins - number of bins for tuning
%
% Outputs:
%   foldEval - NeuronTable structure with evaluation information from fold
%       Each row corresponds to a neuron's evaluation over both workspaces
%       .{model_name}_eval - pseudo-R2 for a given model to that neuron
%   foldTuning - NeuronTable structure with tuning weight information from fold
%       Each row corresponds to a neuron evaluated in one of the workspaces
%       Columns past the header columns correspond to evaluation criteria:
%       .{model_name}_*PD - classical PD for given model's predicted firing rates
%       .{model_name}_*curve - empirical tuning curve for given model's
%           predicted firing rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up
    % default parameters
    model_eval_metric = 'pr2';
    glm_params = {};
    model_names = {};
    num_tuning_bins = 16;
    unit_guide = [];
    crossvalID = [];
    if nargin > 2
        assignParams(who,params);
    end % overwrite parameters

    % check inputs
    assert(~isempty(glm_params),'Must pass in glm parameters')
    assert(~isempty(model_names),'Must pass in model names')
    assert(length(model_names) == length(glm_params) + 1,'Model names must have one more element than glm_params')
    assert(~isempty(crossvalID),'Must pass in crossvalID')

%% Fit models
    % set up parameters for models
    glm_info = cell(1,length(model_names)-1);
    glm_info_within = cell(2,length(model_names)-1);
    for modelnum = 1:length(model_names)-1
        [~,glm_info{modelnum}] = getModel(td_train,glm_params{modelnum});
        
        % get models for training on individual workspaces
        for spacenum = 1:2
            td_train_space = [];
            
            td_train_fields = fieldnames(td_train);
            for i_field = 1:numel(td_train_fields)
                if(length(td_train.(td_train_fields{i_field})) == length(train_spacenum))
                    td_train_space.(td_train_fields{i_field}) = td_train.(td_train_fields{i_field})(train_spacenum==spacenum,:);
                else
                    td_train_space.(td_train_fields{i_field}) = td_train.(td_train_fields{i_field});
                end
            end

            [~,glm_info_within{spacenum,modelnum}] = getModel(td_train_space,glm_params{modelnum});
        end
    end

    % Predict firing rates
    td_test_within = td_test;
    for modelnum = 1:length(model_names)-1
        for spacenum = 1:2
            td_test{spacenum} = getModel(td_test{spacenum},glm_info{modelnum});
            td_test_within{spacenum} = getModel(td_test_within{spacenum},glm_info_within{spacenum,modelnum});
        end
    end

    % Evaluate model fits and add to foldEval table
    foldEval = makeNeuronTableStarter(td_train,struct('out_signal_names',unit_guide,'meta',struct('crossvalID',crossvalID)));
    model_eval = cell(1,length(model_names)-1);
    [space_model_eval,space_model_eval_within] = deal(cell(2,length(model_names)-1));
    eval_params = glm_info;
    for modelnum = 1:length(model_names)-1
        eval_params{modelnum}.eval_metric = model_eval_metric;
        eval_params{modelnum}.num_boots = 1;
        % combine td_test{1} and td_test{2} (and anymore)
        td_test_comb = td_test{1};
        td_field_names = fieldnames(td_test{1});
        for i_td = 2:numel(td_test)
            for i_field = 1:numel(td_field_names)
                if(~ischar(td_test_comb.(td_field_names{i_field})) && length(td_test_comb.(td_field_names{i_field})) > 5)
                    td_test_comb.(td_field_names{i_field}) = [td_test_comb.(td_field_names{i_field}); td_test{i_td}.(td_field_names{i_field})];
                end
            end
        end
        model_eval{modelnum} = array2table(...
            squeeze(evalModel(td_test_comb,eval_params{modelnum}))',...
            'VariableNames',{strcat(model_names{modelnum},'_eval')});
        model_eval{modelnum}.Properties.VariableDescriptions = {'linear'};

        % get evals for individual spaces
        for spacenum = 1:2
            space_model_eval{spacenum,modelnum} = array2table(...
                squeeze(evalModel(td_test{spacenum},eval_params{modelnum}))',...
                'VariableNames',{sprintf('%s_space%d_eval',model_names{modelnum},spacenum)});
            space_model_eval{spacenum,modelnum}.Properties.VariableDescriptions = {'linear'};
            
            space_model_eval_within{spacenum,modelnum} = array2table(...
                squeeze(evalModel(td_test_within{spacenum},eval_params{modelnum}))',...
                'VariableNames',{sprintf('%s_space%d_within_eval',model_names{modelnum},spacenum)});
            space_model_eval_within{spacenum,modelnum}.Properties.VariableDescriptions = {'linear'};
        end
    end
    foldEval = horzcat(foldEval, model_eval{:}, space_model_eval{:}, space_model_eval_within{:});

%% Get extrinsic test tuning (to calculate later quantities from)
    tempTuningTable = cell(2,1);
    for spacenum = 1:2
        for modelnum = 1:length(model_names)
            % get tuning weights for each model
            pdParams = struct(...
                'out_signals',model_names(modelnum),...
                'prefix',model_names{modelnum},...
                'out_signal_names',unit_guide,...
                'bootForTuning',false,...
                'num_boots',50,...
                'verbose',false,...
                'meta',struct('spaceNum',spacenum,'crossvalID',crossvalID));
            temp_pdTable = getTDClassicalPDs(td_test{spacenum},pdParams);
            % temp_pdTable = getTDPDs(td_test{spacenum},pdParams);

            tuningParams = struct(...
                'out_signals',model_names(modelnum),...
                'prefix',model_names{modelnum},...
                'out_signal_names',unit_guide,...
                'calc_CIs',false,...
                'num_bins',num_tuning_bins,...
                'meta',struct('spaceNum',spacenum,'crossvalID',crossvalID));
            temp_tuning_table = getTuningCurves(td_test{spacenum},tuningParams);
            temp_table = join(temp_pdTable,temp_tuning_table);

            % append table to full tuning table for space
            if modelnum == 1
                tempTuningTable{spacenum} = temp_table;
            else
                tempTuningTable{spacenum} = join(tempTuningTable{spacenum}, temp_table);
            end
        end
    end
    % smoosh space tables together
    foldTuning = vertcat(tempTuningTable{:});

end
