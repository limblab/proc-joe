function results = twoPlanesEncoders(td_list,robot_height,params)
%% Split up trial data and preprocess
    % Split td into different workspaces each entry in td_list is a
    % different workspace. Each 2D file has the robot at a different height,
    % specified by robot_height as a value in cm. Larger numbers typically mean higher up (closer to monkey's head than floor).
    
    % also make sure we have balanced workspaces (slightly biases us towards early trials, but this isn't too bad)
    if(robot_height(1) > robot_height(2))
        td_high = td_list{1};
        td_low = td_list{2};
    else
        td_high = td_list{2};
        td_low = td_list{1};
    end
    
    min_dur = min(length(td_high.pos),length(td_low.pos));
    for i=1:2
        switch i
            case 1
                td_temp = td_high;
            case 2
                td_temp = td_low;
        end
        
        td_fields = fieldnames(td_temp);
        pos_length = length(td_temp.pos);
        for i_field = 1:numel(td_fields)
            if(size(td_temp.(td_fields{i_field}),1) == pos_length)
                td_temp.(td_fields{i_field}) = td_temp.(td_fields{i_field})(1:min_dur,:);
            end
        end
        
        switch i
            case 1
                td_high = td_temp;
            case 2
                td_low = td_temp;
        end
        clear td_temp;
    end

%% Set up model variables
    num_folds = 5; % 5 is default number of folds, no need to pass in
    num_repeats = 20; % 20 is default number of repeats, no need to pass in
    num_tuning_bins = 16;
    num_musc_pcs = 5;
    model_type = 'glm';
    model_aliases = {'ext','ego','musc','handelbow'};
    arrayname = 'S1';
    get_tuning_curves = false;
    assignParams(who,params);
    neural_signals = [arrayname '_FR'];
    unit_guide = td_high(1).([arrayname '_unit_guide']);

    model_names = [strcat(model_type,'_',model_aliases,'_model') {neural_signals}];
    num_models = length(model_names);
    
    % set up glm parameters
    glm_params = cell(num_models-1,1);
    for modelnum = 1:num_models-1
        switch model_aliases{modelnum}
        case 'ext'
%             markername = 'Marker_1';
            markername = 'hand2';
            [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_high.dlc_pos_names);
            assert(all(point_exists),'Hand marker does not exist?')
            glm_params{modelnum} = struct(...
                'model_type',model_type,...
                'model_name',[model_aliases{modelnum} '_model'],...
                'in_signals',{{'dlc_pos',marker_hand_idx;'dlc_vel',marker_hand_idx}},...
                'out_signals',neural_signals);
        case 'handelbow'
            % indices for cartesian hand coordinates
            markername = 'hand2';
            [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_high.dlc_pos_names);
            assert(all(point_exists),'Hand marker does not exist?')

            markername = 'elbow1';
            [point_exists,marker_elbow_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_low.dlc_pos_names);
            assert(all(point_exists),'Elbow marker does not exist?')

            glm_params{modelnum} = struct(...
                'model_type',model_type,...
                'model_name',[model_aliases{modelnum} '_model'],...
                'in_signals',{{'dlc_pos',[marker_hand_idx marker_elbow_idx];'dlc_vel',[marker_hand_idx marker_elbow_idx]}},...
                'out_signals',neural_signals);

        otherwise
            error('Unrecognized model_alias')
        end
    end

%% Get comparison of actual tuning curves with various modeled tuning curves
    pdTables = cell(2,num_models);
    tuning_curves = cell(2,num_models);
    td_tuning = [];
    if get_tuning_curves
        % use K-fold crossvalidation to get neural predictions from each model for tuning curves and PDs
        indices = crossvalind('Kfold',length(td_high.pos),num_folds);
        td_test_high = cell(1,num_folds);
        td_test_low = cell(1,num_folds);

        % do the crossval
        for foldctr = 1:num_folds
            % split into testing and training
            test_idx = (indices==foldctr);
            train_idx = ~test_idx;
            td_train = [];

            % split into td_train and td_test, only use in_signals and
            % out_signals. Also grab handle velocity for td_test so that we
            % can compute PDs using handle kinematics
            
            for i_name = 1:size(glm_params{modelnum}.in_signals,2)
                for i_idx = 1:numel(glm_params{modelnum}.in_signals{i_name,2})
                    td_train.(glm_params{modelnum}.in_signals{i_name,1}) = ...
                        [td_high.(glm_params{modelnum}.in_signals{i_name,1})(train_idx,:); td_low.(glm_params{modelnum}.in_signals{i_name,1})(train_idx,:)];
                    td_test_high{1,foldctr}.(glm_params{modelnum}.in_signals{i_name,1}) = td_high.(glm_params{modelnum}.in_signals{i_name,1})(test_idx,:);
                    td_test_low{1,foldctr}.(glm_params{modelnum}.in_signals{i_name,1}) = td_low.(glm_params{modelnum}.in_signals{i_name,1})(test_idx,:);
                    td_test_high{1,foldctr}.vel = td_high.vel(test_idx,:);
                    td_test_low{1,foldctr}.vel = td_low.vel(test_idx,:);
                    
                    % meta info
                    td_test_high{1,foldctr}.monkey = td_high.monkey;
                    td_test_low{1,foldctr}.monkey = td_low.monkey;
                    
                    td_test_high{1,foldctr}.task = td_high.task;
                    td_test_low{1,foldctr}.task = td_low.task;
                    
                    td_test_high{1,foldctr}.date = td_high.date;
                    td_test_low{1,foldctr}.date = td_low.date;
                end
            end
            % get out signals
            td_train.(glm_params{modelnum}.out_signals) = ...
                [td_high.(glm_params{modelnum}.out_signals)(train_idx,:); td_low.(glm_params{modelnum}.out_signals)(train_idx,:)];
            td_test_high{1,foldctr}.(glm_params{modelnum}.out_signals) = td_high.(glm_params{modelnum}.out_signals)(test_idx,:); 
            td_test_low{1,foldctr}.(glm_params{modelnum}.out_signals) = td_low.(glm_params{modelnum}.out_signals)(test_idx,:);
            
            % set bin_size as well
            td_train.bin_size = td_high.bin_size;
            td_test_high{1,foldctr}.bin_size = td_low.bin_size;
            td_test_low{1,foldctr}.bin_size = td_low.bin_size;
            
            % Fit models on training data
            for modelnum = 1:num_models-1
                [~,glm_info] = getModel(td_train,glm_params{modelnum});

                % predict firing rates for td_test
                td_test_high{1,foldctr} = getModel(td_test_high{1,foldctr},glm_info);
                td_test_low{1,foldctr} = getModel(td_test_low{1,foldctr},glm_info);
            end
        end
        
        td_test_high = horzcat(td_test_high{1,:});
        td_test_low = horzcat(td_test_low{1,:});
        
        % get PDs and tuning curves
        for modelnum = 1:num_models
            % First PDs
            pd_params = struct(...
                'out_signals',model_names{modelnum},...
                'out_signal_names',unit_guide,...
                'do_plot',false,...
                'meta',struct('spaceNum',0));
            pdTables{1,modelnum} = getTDClassicalPDs(td_test_high,pd_params);
            pdTables{2,modelnum} = getTDClassicalPDs(td_test_low,pd_params);
            
            
            tuning_params = struct(...
                'out_signals',model_names{modelnum},...
                'out_signal_names',unit_guide,...
                'num_bins',num_tuning_bins,...
                'meta',struct('spaceNum',0));
            tuning_curves{1,modelnum} = getTuningCurves(td_test_high,tuning_params);
            tuning_curves{2,modelnum} = getTuningCurves(td_test_low,tuning_params);
        end
    end

%% Cross-validate models of neural data
    % Get crossval info
    crossval_params = struct(...
        'model_names',{model_names},...
        'glm_params',{glm_params},...
        'num_folds',num_folds,...
        'num_repeats',num_repeats,...
        'crossval_lookup',[],...
        'unit_guide',unit_guide,...
        'num_tuning_bins',num_tuning_bins);
    [crossEval,crossTuning,crossval_lookup] = analyzeTwoPlanes({td_high,td_low},crossval_params);

%% create return struct
    % for cross validation plots
    results.crossEval = crossEval;
    results.crossTuning = crossTuning;

    % get names of neurons
    signalIDs = unit_guide;
    
    if get_tuning_curves
        results.tuning_curves = tuning_curves;
        results.pdTables = pdTables;
        results.isTuned = pdTables{1,end}.velTuned & pdTables{2,end}.velTuned;
        results.tunedNeurons = signalIDs(results.isTuned,:);
        
        % for showing predictive capability
        results.td_tuning = td_tuning;
    end

    % get parameters
    results.params.num_folds = num_folds;
    results.params.num_repeats = num_repeats;
    results.params.num_tuning_bins = num_tuning_bins;
    results.params.model_type = model_type;
    results.params.model_aliases = model_aliases;
    results.params.model_names = model_names;
    results.params.num_musc_pcs = num_musc_pcs;
    results.params.neural_signals = neural_signals;
    results.params.glm_params = glm_params;
    results.params.crossval_lookup = crossval_lookup;
