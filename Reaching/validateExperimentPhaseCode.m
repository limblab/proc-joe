%% Set up meta info and load trial data
    clear; clc;
    if ispc
        folderpath = 'D:\Lab\Data\DLC_videos\Han_20201204_rwFreeReach\neural-data\';
    else
        folderpath = '/data/raeed/project-data/limblab/s1-kinematics';
    end
    
    % load data
    file_info = dir(fullfile(folderpath,'*td*'));
    annotation_file_info = dir(fullfile(folderpath,'*annotation*'));
    
    filenames = horzcat({file_info.name})';
    
    bin_size = 0.05; % s

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
    
%% Loop through trial data files to get masks for cleaning up. Compare to annotated data
   
    filenum = 1;
    % Load data
    load(fullfile(file_info(filenum).folder,file_info(filenum).name));

    keep_mask_all = ones(size(td_list));
    
    for i_td = 2%:numel(td_list) % each entry in td_list is a different trial_data for a different experiment (free reach vs 2D random walk for example)
        % load annotation data
        if(strcmpi(task_list{i_td},'RT3D') == 1)
            load(fullfile(annotation_file_info(1).folder,annotation_file_info(1).name));
        else
            load(fullfile(annotation_file_info(2).folder,annotation_file_info(2).name));
        end
        
        % resample trial data to appropriate bin size
        if(td_list{i_td}.bin_size <= bin_size)
            td_list{i_td} = binTD(td_list{i_td},bin_size/td_list{i_td}.bin_size);
        else
            warning('td bin size is larger than desired bin size');
        end

        if(isfield(td_list{i_td},'dlc_pos'))
            % set origin as shoulder position at t=0
            td_list{i_td} = setOriginAsShoulder(td_list{i_td},0); % use fixed position (t=0) or set shoulder as 0 for each data point.
            % get marker velocity
            td_list{i_td} = getDifferential(td_list{i_td},struct('signals','dlc_pos','alias','dlc_vel'));
            % remove time points where dlc tracking is bad
            dlc_idx = [find((strcmpi(td_list{i_td}.dlc_pos_names,['hand2','_x']))),...
                        find((strcmpi(td_list{i_td}.dlc_pos_names,['hand2','_y']))),...
                        find((strcmpi(td_list{i_td}.dlc_pos_names,['hand2','_z']))),...
                        find((strcmpi(td_list{i_td}.dlc_pos_names,['elbow1','_x']))),...
                        find((strcmpi(td_list{i_td}.dlc_pos_names,['elbow1','_y']))),...
                        find((strcmpi(td_list{i_td}.dlc_pos_names,['elbow1','_y'])))];
            bad_points = any(isnan(td_list{i_td}.dlc_pos(:,:)),2) | any(isnan(td_list{i_td}.dlc_vel(:,:)),2);
            td_names = fieldnames(td_list{i_td});

            fprintf('Removed %.2f%% percent of trials because of missing markers\n',sum(bad_points)/numel(bad_points)*100)
        end
        % remove trials where monkey's arm position is out of the
        % "workspace"
        %td_list{i_td}
        %task_list{i_td}
        [~,experiment_phase_mask,spd_mask,pos_force_mask] = getExperimentPhase(td_list{i_td},task_list{i_td});
        
        
        % plot speed of hand marker
        markername = 'hand2';
        dlc_idx = [find(strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_x'])),find(strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_y'])),...
            find(strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_z']))];
        marker_spd = sqrt(sum(td_list{i_td}.dlc_vel(:,dlc_idx).^2,2));
        figure();
        ax_1=subplot(2,1,1);
        plot(marker_spd);
        
        % plot experiment phase mask
        ax_2=subplot(2,1,2); hold on
        plot(experiment_phase_mask-0.01,'color',getColorFromList(1,1),'linewidth',2)
        % plot bad point mask
        plot(~bad_points + 0.01,'color',getColorFromList(1,0),'linewidth',2);
%         plot annotation labels
        plot(experiment_phase_annotation(end-numel(experiment_phase_mask)+1:end),'k','linewidth',2)
        l=legend('Code','Missing markers','Annotation'); set(l,'box','off');
        linkaxes([ax_1,ax_2],'x');
        
    end
    
    
%% get confusion matrix
    num_points = sum(experiment_phase_mask);
    exp_phase_annotation_align = experiment_phase_annotation(end-numel(experiment_phase_mask)+1:end)';
    
    exp_phase_annotation_align = exp_phase_annotation_align(~bad_points);
    exp_phase_mask_use = experiment_phase_mask(~bad_points);
    
    TN = sum(exp_phase_annotation_align == 0 & exp_phase_mask_use == 0)/sum(exp_phase_annotation_align==0)
    TP = sum(exp_phase_annotation_align == 1 & exp_phase_mask_use == 1)/sum(exp_phase_annotation_align==1)
    FN = sum(exp_phase_annotation_align == 1 & exp_phase_mask_use == 0)/sum(exp_phase_annotation_align==1)
    FP = sum(exp_phase_annotation_align == 0 & exp_phase_mask_use == 1)/sum(exp_phase_annotation_align==0)
    