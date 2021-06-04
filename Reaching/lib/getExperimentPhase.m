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

    % remove data points where speed is too high for hand and elbow
    markernames = {'hand2','elbow1'};
    for marker = markernames
        dlc_idx = [find(strcmpi(trial_data.dlc_pos_names,[marker{1},'_x'])),find(strcmpi(trial_data.dlc_pos_names,[marker{1},'_y'])),...
            find(strcmpi(trial_data.dlc_pos_names,[marker{1},'_z']))];
        marker_spd = sqrt(sum(trial_data.dlc_vel(:,dlc_idx).^2,2));
        max_spd = mean(marker_spd(keep_mask_spd==1))+6*std(marker_spd(keep_mask_spd==1));
        keep_mask_spd = keep_mask_spd & marker_spd < max_spd;
    end
    
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
    % for RW task, check if z-force is non-zero. 
    % for 3D reaching task, check position of hand in workspace
    
    if(strcmpi(task,'RW')==1 || strcmpi(task,'RT') == 1)
        
        z_force = trial_data.force(:,3);
        z_force_thresh = 0.4; % look at histograms of force when the speed is low
        
         % determine whether each point is above the mean marker_spd -
        force_above_thresh = (z_force > z_force_thresh);
        % take moving average to look at segment of data
        mov_mean_force = movmean(force_above_thresh, scan_len);
        % keep segments where enough of data is above mean_speed
        keep_mask_pos = mov_mean_force > prop_above_thresh;
        
        % remove data points where z-value is too far from median
        markername = 'hand2';
        dlc_idx = find(strcmpi(trial_data.dlc_pos_names,[markername,'_z']));
        dlc_pos = trial_data.dlc_pos(:,dlc_idx);
        keep_mask_pos = keep_mask_pos & abs(dlc_pos - median(dlc_pos)) < 6;
        
    elseif(strcmpi(task,'RT3D')==1)
        % remove if hand position is outside predefined cube.
        if(strcmpi(trial_data.monkey,'Han'))
            x_upper_bound = 35; % remove if larger than this
            x_lower_bound = -35; %remove if smaller than this
            y_upper_bound = 30; % remove if smaller than this
            y_lower_bound = -10;
            z_upper_bound = 50;
            z_lower_bound = -20; % remove if larger than this
        elseif(strcmpi(trial_data.monkey,'Crackle'))       
            x_upper_bound = 35; % remove if larger than this, and others
            x_lower_bound = -35; %remove if smaller than this
            y_upper_bound = 30; % remove if smaller than this, and others
            y_lower_bound = -17;
            z_upper_bound = 35;
            z_lower_bound = -15; % remove if larger than this, and others
        else
           warning('position bounds for 3D task is not implemented for this monkey.');
           x_upper_bound = Inf;
           x_lower_bound = -Inf;
           y_upper_bound = Inf;
           y_lower_bound = -Inf;
           z_upper_bound = Inf;
           z_lower_bound = -Inf;
        end
        
         markername = 'hand2';
        dlc_idx = [find(strcmpi(trial_data.dlc_pos_names,[markername,'_x'])),find(strcmpi(trial_data.dlc_pos_names,[markername,'_y'])),...
            find(strcmpi(trial_data.dlc_pos_names,[markername,'_z']))];
        dlc_pos = trial_data.dlc_pos(:,dlc_idx);
        dlc_pos = dlc_pos - median(dlc_pos,1);
        
        keep_mask_pos = dlc_pos(:,1) <= x_upper_bound & dlc_pos(:,1) >= x_lower_bound &...
                dlc_pos(:,2) <= y_upper_bound & dlc_pos(:,2) >= y_lower_bound & ...
                dlc_pos(:,3) <= z_upper_bound & dlc_pos(:,3) >= z_lower_bound;
        
        % remove if ELBOW position is outside predefined cube.
        if(strcmpi(trial_data.monkey,'Han'))
            x_upper_bound = Inf; % remove if larger than this
            x_lower_bound = -Inf; %remove if smaller than this
            y_upper_bound = Inf; % remove if smaller than this
            y_lower_bound = -Inf;
            z_upper_bound = Inf;
            z_lower_bound = -Inf; % remove if larger than this
        elseif(strcmpi(trial_data.monkey,'Crackle'))       
            x_upper_bound = 35; % remove if larger than this, and others
            x_lower_bound = -20; %remove if smaller than this
            y_upper_bound = Inf; % remove if smaller than this, and others
            y_lower_bound = -Inf;
            z_upper_bound = Inf;
            z_lower_bound = -Inf; % remove if larger than this, and others
        else
           warning('position bounds for 3D task is not implemented for this monkey.');
           x_upper_bound = Inf;
           x_lower_bound = -Inf;
           y_upper_bound = Inf;
           y_lower_bound = -Inf;
           z_upper_bound = Inf;
           z_lower_bound = -Inf;
        end
        
        markername = 'elbow1';
        dlc_idx = [find(strcmpi(trial_data.dlc_pos_names,[markername,'_x'])),find(strcmpi(trial_data.dlc_pos_names,[markername,'_y'])),...
            find(strcmpi(trial_data.dlc_pos_names,[markername,'_z']))];
        dlc_pos = trial_data.dlc_pos(:,dlc_idx);
        dlc_pos = dlc_pos - median(dlc_pos,1);
        
        keep_mask_pos = keep_mask_pos & dlc_pos(:,1) <= x_upper_bound & dlc_pos(:,1) >= x_lower_bound &...
                dlc_pos(:,2) <= y_upper_bound & dlc_pos(:,2) >= y_lower_bound & ...
                dlc_pos(:,3) <= z_upper_bound & dlc_pos(:,3) >= z_lower_bound;   
        % only keep positions within 99% bound for each
        % direction....Removing outliers is important,
            
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
