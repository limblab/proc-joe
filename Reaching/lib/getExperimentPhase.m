function [trial_data,f] = getExperimentPhase(trial_data,task)

    % finds time points when monkey is doing the task and returns a trial
    % data containing experiment data only. Also adds field to trial_data
    % that contains which idx (effectively timestamp) corresponds to each
    % data point
    
    % for both tasks, remove spans of time when hand speed is too low (representing when monkey is not working)
    
    markername = 'hand2';
    dlc_idx = [find(strcmpi(trial_data.dlc_pos_names,[markername,'_x'])),find(strcmpi(trial_data.dlc_pos_names,[markername,'_y'])),...
        find(strcmpi(trial_data.dlc_pos_names,[markername,'_z']))];
    dlc_pos = trial_data.dlc_pos(:,dlc_idx);
    
    % Find sections of data (scan_len long) where a large enough proportion
    % (prop above thresh) of speed data is above the mean speed. 
    % 
    marker_spd = sqrt(sum(trial_data.dlc_vel(:,dlc_idx).^2,2));
        exp_phase_spd_thresh = mean(marker_spd);
    scan_len = 40; %frames/bins
    prop_above_thresh = 0.15; %If 60% of the points in this section were considered as 1

    % determine whether each point is above the mean marker_spd -
    spd_above_thresh = (marker_spd > exp_phase_spd_thresh);
    % take moving average to look at segment of data
    mov_mean_spd_above_thresh = movmean(spd_above_thresh, scan_len);
    % keep segments where enough of data is above mean_speed
    keep_mask_spd = mov_mean_spd_above_thresh > prop_above_thresh;

    % remove data points where speed is too high
    max_spd = mean(marker_spd(keep_mask_spd==1))+6*std(marker_spd(keep_mask_spd==1));
    keep_mask_spd = keep_mask_spd & marker_spd < max_spd;

    %Plot the actual speed and whether we kept the data point (1=yes,0=no)
    x_data = 0:trial_data.bin_size:(length(marker_spd))*trial_data.bin_size;
    f=figure;
    f.Name = [trial_data.monkey, '_', trial_data.date, '_', task, '_speedFilter'];
    plot(x_data(1:length(marker_spd)),marker_spd);
    hold on
    plot(x_data(1:length(keep_mask_spd)),keep_mask_spd*mean(marker_spd),'LineWidth',1);
    hold off 
    xlim([min(x_data), max(x_data)]);
    l=legend('Actual Marker Speed','Predicted Mask');
    set(l,'box','off');
    set(gca,'Fontsize',14);
    formatForLee(f);
    xlabel('Time (s)');
    ylabel('Speed (cm/s)');
    
    % remove data points if outside of position thresholds
    % for RW task, check if monkey's hand is reasonably close to the handle
    % (use hand1,2,or 3 marker). Also make sure handle is not still for
    % extended periods of time
    % for 3D reaching task, check position of hand in workspace
    
    if(strcmpi(task,'RW')==1 || strcmpi(task,'RT') == 1)
        % get transformation from DLC coordinates to robot coordinates
        % (rotation and offset)
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
        keep_mask_pos = sqrt(sum((dlc_rot(:,1:2)-handle_pos).^2,2)) < 5; % 5 was chosen again by looking at a histogram of distances for an example dataset
        
    elseif(strcmpi(task,'RT3D')==1)
        % remove if hand position is outside predefined cube.
        if(strcmpi(trial_data.monkey,'Han'))
            x_upper_bound = 16; % remove if larger than this
            x_lower_bound = -40; %remove if smaller than this
            y_upper_bound = 40; % remove if smaller than this
            y_lower_bound = 15;
            z_upper_bound = Inf;
            z_lower_bound = -10; % remove if larger than this
        elseif(strcmpi(trial_data.monkey,'Crackle'))       
            x_upper_bound = 10; % remove if larger than this, and others
            x_lower_bound = -30; %remove if smaller than this
            y_upper_bound = 45; % remove if smaller than this, and others
            y_lower_bound = 20;
            z_upper_bound = 40;
            z_lower_bound = -10; % remove if larger than this, and others
        else
           warning('position bounds for 3D task is not implemented for this monkey.');
           x_upper_bound = Inf;
           x_lower_bound = -Inf;
           y_upper_bound = Inf;
           y_lower_bound = -Inf;
           z_upper_bound = Inf;
           z_lower_bound = -Inf;
        end
        
        keep_mask_pos = trial_data.dlc_pos(:,dlc_idx(1)) <= x_upper_bound & trial_data.dlc_pos(:,dlc_idx(1)) >= x_lower_bound &...
                trial_data.dlc_pos(:,dlc_idx(2)) <= y_upper_bound & trial_data.dlc_pos(:,dlc_idx(2)) >= y_lower_bound & ...
                trial_data.dlc_pos(:,dlc_idx(3)) <= z_upper_bound & trial_data.dlc_pos(:,dlc_idx(3)) >= z_lower_bound;
        
    end %end of either process the RT2D task or the RT3D task to experiment-only data
    
    
    % remove points
    td_fields = fieldnames(trial_data);
    pos_len = length(trial_data.dlc_pos);
    for i_field = 1:numel(td_fields)
        if(length(trial_data.(td_fields{i_field})) == pos_len)
            trial_data.(td_fields{i_field}) = trial_data.(td_fields{i_field})(keep_mask_spd==1 & keep_mask_pos == 1,:);
        end
    end
    trial_data.trial_idx = find(keep_mask_spd & keep_mask_pos);
    
    
end
