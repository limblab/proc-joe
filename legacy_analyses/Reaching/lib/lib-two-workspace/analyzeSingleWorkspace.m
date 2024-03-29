%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [crossEval, crossTuning, crossvalLookup] = analyzeSingleWorkspace(trial_data,params)
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
function [crossEval, crossTuning, crossvalLookup] = analyzeSingleWorkspace(td,params)
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

    % loop over num repeats
    if verbose
        repeat_timer = tic;
        fprintf('Starting %dx%d-fold crossvalidation at time %f\n',num_repeats,num_folds,toc(repeat_timer));
    end
    for repeatctr = 1:num_repeats
        % get fold indices
        indices = crossvalind('Kfold',length(td.pos),num_folds);

        % initialize temporary fold evaluation structure
        [foldEval,foldTuning,foldCrossvalLookup] = deal(cell(num_folds,1));

        % loop over number of folds
        if verbose
            fold_timer = tic;
        end
        for foldctr = 1:num_folds
            % Get test and training indices for this fold
            % check if there's a crossval lookup table
            td_test = [];
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
                                td.(glm_params{modelnum}.in_signals{i_name,1})(train_idx,:);
                            td_test.(glm_params{modelnum}.in_signals{i_name,1}) = td.(glm_params{modelnum}.in_signals{i_name,1})(test_idx,:);
                            % meta info....
                            td_train.monkey = td.monkey;
                            td_test.monkey = td.monkey;
            
                            td_train.task = td.task;
                            td_test.task = td.task;

                            td_train.date = td.date;
                            td_test.date = td.date;
                        end
                    end
                end
                
                % get out signals
                td_train.(glm_params{1}.out_signals) = ...
                    td.(glm_params{1}.out_signals)(train_idx,:);
                td_test.(glm_params{1}.out_signals) = td.(glm_params{1}.out_signals)(test_idx,:); 

                % set bin_size as well
                td_train.bin_size = td.bin_size;
                td_test.bin_size = td.bin_size;
                
                % set spacenum (1 only for a single workspace)
                train_spacenum = ones(sum(train_idx),1);
            else
                error('crossval lookup is not empty. This is not implemented');
            end

            % analyze fold to get model evaluations
            if exist('params','var')
                params.crossvalID = uint16([repeatctr foldctr]);
            else
                params = struct('crossvalID',uint16([repeatctr foldctr]));
            end
            [foldEval{foldctr},foldTuning{foldctr}] = analyzeFold(td_train,td_test,train_spacenum,params);

            % get test trialIDs
            trialID = cat(1,find(indices==foldctr)',find(indices==foldctr)');
            spaceNum = zeros(size(trialID)) + [1]; % 1 only for single workspace
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
    
    glm_info = cell(1,length(model_names)-1);
    glm_info_within = cell(2,length(model_names)-1);
    for modelnum = 1:length(model_names)-1
        [~,glm_info{modelnum}] = getModel(td_train,glm_params{modelnum});
    end
    
    % Predict firing rates
    for modelnum = 1:length(model_names)-1
        td_test = getModel(td_test,glm_info{modelnum});
    end

    % Evaluate model fits and add to foldEval table
    foldEval = makeNeuronTableStarter(td_train,struct('out_signal_names',unit_guide,'meta',struct('crossvalID',crossvalID)));
    model_eval = cell(1,length(model_names)-1);
    eval_params = glm_info;
    for modelnum = 1:length(model_names)-1
        eval_params{modelnum}.eval_metric = model_eval_metric;
        eval_params{modelnum}.num_boots = 1;
        
        model_eval{modelnum} = array2table(...
            squeeze(evalModel(td_test,eval_params{modelnum}))',...
            'VariableNames',{strcat(model_names{modelnum},'_eval')});
        model_eval{modelnum}.Properties.VariableDescriptions = {'linear'};

    end
    foldEval = horzcat(foldEval, model_eval{:});

%% Get extrinsic test tuning (to calculate later quantities from)
    foldTuning = {};

end
