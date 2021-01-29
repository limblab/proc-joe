function [trial_data] = getExperimentPhase(trial_data,task)

    % finds time points when monkey is doing the task and returns a trial
    % data containing experiment data only. Also adds field to trial_data
    % that contains which idx (effectively timestamp) corresponds to each
    % data point
    
    
    % for RW task, check if monkey's hand is reasonably close to the handle
    % (use hand1,2,or 3 marker). Also make sure handle is not still for
    % extended periods of time
    
    % for 3D reaching task... (not yet implemented
    
    if(strcmpi(task,'RW')==1 || strcmpi(task,'RT') == 1)
        % get transformation from DLC coordinates to robot coordinates
        % (rotation and offset)
        markername = 'hand3';
        dlc_idx = [find(strcmpi(trial_data.dlc_pos_names,[markername,'_x'])),find(strcmpi(trial_data.dlc_pos_names,[markername,'_y'])),...
            find(strcmpi(trial_data.dlc_pos_names,[markername,'_z']))];
        dlc_pos = trial_data.dlc_pos(:,dlc_idx);
        handle_pos = [trial_data.pos]; 
        
        % use only times when handle is moving for fit
        handle_mask = sqrt(sum(trial_data.vel.^2,2)) > 1; % 1 was chosen by looking at a histogram of speeds
        handle_pos = handle_pos(handle_mask==1,:);
        dlc_pos = dlc_pos(handle_mask==1,:);
        % subtrack mean
        mean_dlc_pos = mean(dlc_pos);
        mean_handle_pos = mean(handle_pos);
        dlc_pos = dlc_pos - mean_dlc_pos;
        handle_pos = handle_pos - mean_handle_pos;
        
        % find rotation matrix
        rot_mat = dlc_pos(:,1:2)\handle_pos;
        rot_mat = rot_mat./sqrt(sum(rot_mat.^2));
        
        % rotate original data
        dlc_rot = trial_data.dlc_pos(:,dlc_idx);
        dlc_rot = dlc_rot - mean_dlc_pos;
        dlc_rot(:,1:2) = dlc_rot(:,1:2)*rot_mat;
        
        handle_pos = [trial_data.pos] - mean_handle_pos; 
        % check to see if hand is within distance bound to handle
        keep_mask = sqrt(sum((dlc_rot(:,1:2)-handle_pos).^2,2)) < 3.5; % 3.5 was chosen again by looking at a histogram of distances for an example dataset
        
        % remove points
        td_fields = fieldnames(trial_data);
        pos_len = length(trial_data.pos);
        for i_field = 1:numel(td_fields)
            if(length(trial_data.(td_fields{i_field})) == pos_len)
                trial_data.(td_fields{i_field}) = trial_data.(td_fields{i_field})(keep_mask==1,:);
            end
        end
        trial_data.trial_idx = find(keep_mask);
        
    elseif(strcmpi(task,'RT3D')==1)
        warning('getExperimentPhase is not implemented for this task. Nothing was done');
    end
    
    
    




end