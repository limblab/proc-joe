%% Set up meta info and load trial data
    clear; clc;
    if ispc
        folderpath = 'D:\Lab\Data\FreeReaching\Crackle_20201203_rwFreeReach\neural-data\';
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
    
    bin_size = 0.05;
    arrayname = 'LeftS1';
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
        end
        td_all{filenum} = td_list;
        task_list_all{filenum} = task_list;
        robot_height_all{filenum} = robot_height;
        clear td_list;
    end
    td_list = td_all{1};

%% for each trial data, get spindle model PDs

    n_spindles = 2000; % number of spindles for each muscle. spindles are proportioned based on getMuscleSpindleProportions

    pd_table = {};
    for i_td = 1:numel(td_list)
        % get muscle names
        opensim_idx = find(~cellfun(@isempty,strfind(td_list{i_td}.opensim_names,'_len')));

        % get mean muscle length and proportion of spindles in each muscle
        prop = getMuscleSpindleProportions(mean(td_list{i_td}.opensim(:,opensim_idx),1,'omitnan'),td_list{i_td}.opensim_names(opensim_idx));
        % get muscle idx for each spindle
        spindle_musc_idx = datasample(1:1:numel(prop),n_spindles,'Replace',true,'Weights',prop);

        % get muscle length velocity 
        opensim_idx = find(~cellfun(@isempty,strfind(td_list{i_td}.opensim_names,'muscVel')));
        musc_vel = td_list{i_td}.opensim(:,opensim_idx);
        musc_vel(musc_vel < 0) = 0;

        % pass musc_vel through power law with coef 0.5
        musc_fr = musc_vel.^0.5;

        % set 90th percentile firing rate to 50 Hz
        musc_fr_90 = prctile(musc_fr,90);
        musc_fr = 50*musc_fr./musc_fr_90;

        % sample from a poisson distribution
        spindle_lambdas = zeros(size(musc_fr,1),sum(n_spindles));

        for i_spindle = 1:sum(n_spindles)
            spindle_lambdas(:,i_spindle) = musc_fr(:,spindle_musc_idx(i_spindle));        
        end

        spindle_fr = poissrnd(spindle_lambdas);

        % store spindle info in td

        td_list{i_td}.spindle_fr = spindle_fr;
        td_list{i_td}.spindle_musc_idx = spindle_musc_idx;
        td_list{i_td}.musc_names = td_list{i_td}.opensim_names(opensim_idx);

        % use linear model to determine PDs for the simulated spindles based on
        % hand vel (as we do for neurons)
        % get hand vel from dlc
        markername = 'hand2';
        dlc_idx = [find(strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_x'])),...
            find(strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_y'])),...
            find(strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_z']))];

        td_list{i_td}.dlc_vel_hand = td_list{i_td}.dlc_vel(:,dlc_idx);
        pd_params = struct(...
                        'out_signals','spindle_fr',...
                        'in_signals','dlc_vel_hand',...
                        'bootForTuning',false,...
                        'num_boots',50,...
                        'verbose',false,...
                        'doing_3D_pd',true);

        pd_table{i_td} = getTDPDs3D(td_list{i_td},pd_params);

        figure();
        rose(pd_table{i_td}.dlc_vel_handPD)

        figure();
        histogram(pd_table{i_td}.dlc_vel_handPD_zAng);
    end
    task_3d_idx = find(strcmpi('RT3D',task_list));
    task_2d_idx = find(strcmpi('RT',task_list));
    
    tuning_data = [];
    tuning_data.velPD_2D = pd_table{task_2d_idx}.dlc_vel_handPD;
    tuning_data.velPD_zAng_2D = pd_table{task_2d_idx}.dlc_vel_handPD_zAng;
    tuning_data.velPD_3D = pd_table{task_3d_idx}.dlc_vel_handPD;
    tuning_data.velPD_zAng_3D = pd_table{task_3d_idx}.dlc_vel_handPD_zAng;
    
%% get distribution of hand velocities (in plane and z-ang) for comparison
    
    hand_vel_planar = {};
    hand_vel_z_ang = {};
    for i_td = 1:numel(td_list)
        markername = 'hand2';
        dlc_idx = [find((strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_x']))),...
                    find((strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_y']))),...
                    find((strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_z'])))];

        hand_vels = td_list{i_td}.dlc_vel(:,dlc_idx);
        
        % get 2D projection
        hand_vel_planar{i_td} = atan2(hand_vels(:,2),hand_vels(:,1));
        
        % get z-ang
        hand_vel_z_ang{i_td} = atan2(hand_vels(:,3),sqrt(sum(hand_vels(:,1:2).^2,2))); % these two methods are equivalent, except the sign...
                    
        figure();
        histogram(hand_vel_planar{i_td});
        figure();
        histogram(hand_vel_z_ang{i_td});
        
    end
    task_3d_idx = find(strcmpi('RT3D',task_list));
    task_2d_idx = find(strcmpi('RT',task_list));
    
    tuning_data = [];
    tuning_data.velPD_2D = hand_vel_planar{task_2d_idx};
    tuning_data.velPD_zAng_2D = hand_vel_z_ang{task_2d_idx};
    tuning_data.velPD_3D = hand_vel_planar{task_3d_idx};
    tuning_data.velPD_zAng_3D = hand_vel_z_ang{task_3d_idx};
