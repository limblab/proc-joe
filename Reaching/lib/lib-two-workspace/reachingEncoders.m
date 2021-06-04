function results = reachingEncoders(td_list,task_list,robot_height,params)
%% Split up trial data and preprocess
    % Split td into different workspaces each entry in td_list is a
    % different workspace. Each 2D file has the robot at a different height,
    % specified by robot_height as a value in cm. Larger numbers typically mean higher up (closer to monkey's head than floor).
    
    % also make sure we have balanced workspaces (slightly biases us towards early trials, but this isn't too bad)
    if(strcmpi(task_list{1},'RT')==1)
        td_plane = td_list{1};
        td_freereach = td_list{2};
    else
        td_plane = td_list{2};
        td_freereach = td_list{1};
    end
    min_dur = Inf;
    for i_td = 1:numel(td_list)
        min_dur = min(min_dur,length(td_list{i_td}.dlc_pos));
    end
    
    for i=1:2
        switch i
            case 1
                td_temp = td_plane;
            case 2
                td_temp = td_freereach;
        end
        
        td_fields = fieldnames(td_temp);
        pos_length = length(td_temp.dlc_pos);
        for i_field = 1:numel(td_fields)
            if(size(td_temp.(td_fields{i_field}),1) == pos_length)
                td_temp.(td_fields{i_field}) = td_temp.(td_fields{i_field})(1:min_dur,:);
            end
        end
        
        switch i
            case 1
                td_plane = td_temp;
            case 2
                td_freereach = td_temp;
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
    unit_guide = td_plane(1).([arrayname '_unit_guide']);

    model_names = [strcat(model_type,'_',model_aliases,'_model') {neural_signals}];
    num_models = length(model_names);
    
    % set up glm parameters
    glm_params = cell(num_models-1,1);
    for modelnum = 1:num_models-1
        switch model_aliases{modelnum}
            case 'ext'
    %             markername = 'Marker_1';
                markername = 'hand2';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_plane.dlc_pos_names);
                assert(all(point_exists),'Hand marker does not exist?')
                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'dlc_pos',marker_hand_idx;'dlc_vel',marker_hand_idx}},...
                    'out_signals',neural_signals);
                
            case 'ext_2D'
    %             markername = 'Marker_1';
                markername = 'hand2';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y'}),td_plane.dlc_pos_names);
                assert(all(point_exists),'Hand marker does not exist?')
                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'dlc_pos',marker_hand_idx;'dlc_vel',marker_hand_idx}},...
                    'out_signals',neural_signals);    
                
            case 'handelbow'
                % indices for cartesian hand coordinates
                markername = 'hand2';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_plane.dlc_pos_names);
                assert(all(point_exists),'Hand marker does not exist?')

                markername = 'elbow1';
                [point_exists,marker_elbow_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_plane.dlc_pos_names);
                assert(all(point_exists),'Elbow marker does not exist?')

                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'dlc_pos',[marker_hand_idx marker_elbow_idx];'dlc_vel',[marker_hand_idx marker_elbow_idx]}},...
                    'out_signals',neural_signals);
                
            case 'handelbow_2D'
                % indices for cartesian hand coordinates
                markername = 'hand2';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y'}),td_plane.dlc_pos_names);
                assert(all(point_exists),'Hand marker does not exist?')

                markername = 'elbow1';
                [point_exists,marker_elbow_idx] = ismember(strcat(markername,'_',{'x','y'}),td_plane.dlc_pos_names);
                assert(all(point_exists),'Elbow marker does not exist?')

                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'dlc_pos',[marker_hand_idx marker_elbow_idx];'dlc_vel',[marker_hand_idx marker_elbow_idx]}},...
                    'out_signals',neural_signals);

            case 'opensim_ext'
                % indices for cartesian hand coordinates
                markername = 'handPos';
                [point_exists,hand_pos_idx] = ismember(strcat({'X','Y','Z'},'_',markername),td_plane.opensim_names);
                assert(all(point_exists),'Hand marker does not exist?')

                markername = 'handVel';
                [point_exists,hand_vel_idx] = ismember(strcat({'X','Y','Z'},'_',markername),td_plane.opensim_names);
                assert(all(point_exists),'Hand marker does not exist?')

                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'opensim',[hand_pos_idx hand_vel_idx]}},...
                    'out_signals',neural_signals);

            case 'opensim_handelbow'
                % indices for cartesian hand coordinates
                markername = 'handPos';
                [point_exists,hand_pos_idx] = ismember(strcat({'X','Y','Z'},'_',markername),td_plane.opensim_names);
                assert(all(point_exists),'Hand marker does not exist?')

                markername = 'handVel';
                [point_exists,hand_vel_idx] = ismember(strcat({'X','Y','Z'},'_',markername),td_plane.opensim_names);
                assert(all(point_exists),'Hand marker does not exist?')

                markername = 'elbowPos';
                [point_exists,elbow_pos_idx] = ismember(strcat({'X','Y','Z'},'_',markername),td_plane.opensim_names);
                assert(all(point_exists),'elbow marker does not exist?')

                markername = 'elbowVel';
                [point_exists,elbow_vel_idx] = ismember(strcat({'X','Y','Z'},'_',markername),td_plane.opensim_names);
                assert(all(point_exists),'elbow marker does not exist?')

                glm_params{modelnum} = struct('model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'opensim',[hand_pos_idx hand_vel_idx elbow_pos_idx elbow_vel_idx]}},...
                    'out_signals',neural_signals);


            case 'joint'
                markername = 'ang';
                [ang_idx] = find(~cellfun(@isempty,strfind(td_plane.opensim_names,strcat('_',markername))));
                if(isempty(ang_idx))
                    assert('joint angles do not exist?');
                end

                markername = 'vel';
                [vel_idx] = find(~cellfun(@isempty,strfind(td_plane.opensim_names,strcat('_',markername))));
                if(isempty(vel_idx))
                    assert('joint angles do not exist?');
                end

                glm_params{modelnum} = struct('model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'opensim',[ang_idx, vel_idx]}},...
                    'out_signals',neural_signals);



            case 'musc'
                % Do PCA on muscle space
                % do PCA on muscles, training on only the training set
                % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
                % PCA is done outside of this function for some reason....hmm

                glm_params{modelnum} = struct('model_type',model_type,...
                                        'model_name',[model_aliases{modelnum} '_model'],...
                                        'in_signals',{{'musc_len_pca',1:num_musc_pcs;'musc_vel_pca',1:num_musc_pcs}},...
                                        'out_signals',neural_signals);
        otherwise
            error('Unrecognized model_alias')
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
        'num_tuning_bins',num_tuning_bins,...
        'get_tuning_curves',get_tuning_curves);
    [crossEval,crossTuning,crossval_lookup,crossval_glm_info] = analyze3DReachVs2DReach({td_plane,td_freereach},crossval_params);

%% create return struct
    % for cross validation plots
    results.crossEval = crossEval;
    results.crossTuning = crossTuning;

    % get names of neurons
    signalIDs = unit_guide;
    

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
    results.glm_info = crossval_glm_info;
