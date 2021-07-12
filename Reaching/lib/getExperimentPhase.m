function [trial_data,keep_mask_all,keep_mask_spd,keep_mask_pos] = getExperimentPhase(trial_data,task)

    % finds time points when monkey is doing the task and returns a trial
    % data containing experiment data only. Also adds field to trial_data
    % that contains which idx (effectively timestamp) corresponds to each
    % data point
    
    % for both tasks, remove spans of time when hand speed is too low (representing when monkey is not working)
    
    markername = 'hand2';
    dlc_idx = [find(strcmpi(trial_data.dlc_pos_names,[markername,'_x'])),find(strcmpi(trial_data.dlc_pos_names,[markername,'_y'])),...
        find(strcmpi(trial_data.dlc_pos_names,[markername,'_z']))];

    
    
    % go through data set and find long periods of low speed, mark as
    % non-experiment
    marker_spd = sqrt(sum(trial_data.dlc_vel(:,dlc_idx).^2,2));

    scan_len = 60; % bins
    marker_spd_app = [zeros(scan_len,1);marker_spd;zeros(scan_len,1)]; % append zeros for rolling average
    marker_spd_app = fillmissing(marker_spd_app,'constant',0);
    
    marker_spd_roll_avg = movmean(marker_spd_app,scan_len,'omitnan');
    
    exp_spd_thresh = 0.5*mean(marker_spd_roll_avg,'omitnan');
    non_exp_spd_thresh = 0.25*mean(marker_spd_roll_avg,'omitnan');
    
    keep_mask_spd = ones(numel(marker_spd_app),1);
    is_exp_flag = 0; % start as non-experiment phase
    for i = 1:numel(marker_spd_roll_avg)
        if(is_exp_flag == 0 && marker_spd_roll_avg(i) > exp_spd_thresh)
            is_exp_flag = 1;
        elseif(is_exp_flag == 1 && marker_spd_roll_avg(i) < non_exp_spd_thresh)
            is_exp_flag = 0;
            keep_mask_spd(i-ceil(scan_len/3):i) = is_exp_flag; % there's a bit of a delay going down
        end
        keep_mask_spd(i) = is_exp_flag;
    end
    % remove appended zeros
    keep_mask_spd = keep_mask_spd(scan_len+1:end-scan_len);
    
    
%     % Find sections of data (scan_len long) where a large enough proportion
%     % (prop above thresh) of speed data is above the mean speed. 
%     % 
%     marker_spd = sqrt(sum(trial_data.dlc_vel(:,dlc_idx).^2,2));
%     exp_phase_spd_thresh = 0.75*mean(marker_spd,'omitnan');
%     scan_len = 40; %bins
%     prop_above_thresh = 0.15; 
% 
%     % determine whether each point is above the mean marker_spd -
%     spd_above_thresh = (marker_spd > exp_phase_spd_thresh);
%     % take moving average to look at segment of data
%     mov_mean_spd_above_thresh = movmean(spd_above_thresh, scan_len,'omitnan');
%     % keep segments where enough of data is above mean_speed
%     keep_mask_spd = mov_mean_spd_above_thresh > prop_above_thresh;

    
    
    
    
    % remove data points where speed is too high for hand and elbow
    markernames = {'hand2','elbow1'};
    for marker = markernames
        dlc_idx = [find(strcmpi(trial_data.dlc_pos_names,[marker{1},'_x'])),find(strcmpi(trial_data.dlc_pos_names,[marker{1},'_y'])),...
            find(strcmpi(trial_data.dlc_pos_names,[marker{1},'_z']))];
        marker_spd = sqrt(sum(trial_data.dlc_vel(:,dlc_idx).^2,2));
        max_spd = mean(marker_spd(keep_mask_spd==1),'omitnan')+8*std(marker_spd(keep_mask_spd==1),'omitnan');
        keep_mask_spd = keep_mask_spd & marker_spd < max_spd;
    end
    
    %Plot the actual speed and whether we kept the data point (1=yes,0=no)
%     x_data = 0:trial_data.bin_size:(length(marker_spd))*trial_data.bin_size;
%     f=figure;
%     f.Name = [trial_data.monkey, '_', trial_data.date, '_', task, '_speedFilter'];
%     plot(x_data(1:length(marker_spd)),marker_spd);
%     hold on
%     plot(x_data(1:length(keep_mask_spd)),keep_mask_spd*mean(marker_spd),'LineWidth',1);
%     hold off 
%     xlim([min(x_data), max(x_data)]);
%     l=legend('Actual Marker Speed','Predicted Mask');
%     set(l,'box','off');
%     set(gca,'Fontsize',14);
%     formatForLee(f);
%     xlabel('Time (s)');
%     ylabel('Speed (cm/s)');
    
    % remove data points if outside of position thresholds
    % for RW task, check if z-force is non-zero. 
    % for 3D reaching task, check position of hand in workspace
    
    if(strcmpi(task,'RW')==1 || strcmpi(task,'RT') == 1)
        % go through data set and find long periods of low speed, mark as
        % non-experiment
        z_force = trial_data.force(:,3);
        exp_z_force_thresh = 0.4; % look at histograms of force when the speed is low
        non_z_force_thresh = 0.2;
        
        scan_len = 60; % bins
        z_force_app = [zeros(scan_len,1);z_force;zeros(scan_len,1)]; % append zeros for rolling average
        z_force_roll_avg = movmean(z_force_app,scan_len,'omitnan');

        exp_spd_thresh = 0.75*mean(z_force_roll_avg,'omitnan');
        non_exp_spd_thresh = 0.25*mean(z_force_roll_avg,'omitnan');

        keep_mask_pos = ones(numel(z_force_app),1);
        is_exp_flag = 0; % start as non-experiment phase
        for i = 1:numel(z_force_roll_avg)
            if(is_exp_flag == 0 && z_force_roll_avg(i) > exp_z_force_thresh)
                is_exp_flag = 1;
            elseif(is_exp_flag == 1 && z_force_roll_avg(i) < non_z_force_thresh)
                is_exp_flag = 0;
                keep_mask_pos(i-ceil(scan_len/3):i) = is_exp_flag; % there's a bit of a delay going down
            end
            keep_mask_pos(i) = is_exp_flag;
        end
        % remove appended zeros
        keep_mask_pos = keep_mask_pos(scan_len+1:end-scan_len);
        
        % if hand position is out of range, remove points
        markername = 'hand2';
        dlc_idx = [find(strcmpi(trial_data.dlc_pos_names,[markername,'_x'])),find(strcmpi(trial_data.dlc_pos_names,[markername,'_y'])),...
            find(strcmpi(trial_data.dlc_pos_names,[markername,'_z']))];
        dlc_pos = trial_data.dlc_pos(:,dlc_idx);
        dlc_pos = dlc_pos - median(dlc_pos,1,'omitnan');
        
        keep_mask_pos = keep_mask_pos & abs(dlc_pos(:,3)) < 7.5; 
        
        
        
    elseif(strcmpi(task,'RT3D')==1 || strcmpi(task,'freeReach')==1)
        if(isfield(trial_data,'ana_var') && 1==1)
            markername = 'elbow1';
            dlc_idx = [find(strcmpi(trial_data.dlc_pos_names,[markername,'_x'])),find(strcmpi(trial_data.dlc_pos_names,[markername,'_y'])),...
                find(strcmpi(trial_data.dlc_pos_names,[markername,'_z']))];
            
            keep_mask_pos = trial_data.ana_var > 500 & trial_data.dlc_pos(:,dlc_idx(1)) >-55 & trial_data.dlc_pos(:,dlc_idx(1)) < 0;
            
            
            
        else % position based...
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
            dlc_pos = dlc_pos - median(dlc_pos,1,'omitnan');

            keep_mask_pos = dlc_pos(:,1) <= x_upper_bound & dlc_pos(:,1) >= x_lower_bound &...
                    dlc_pos(:,2) <= y_upper_bound & dlc_pos(:,2) >= y_lower_bound & ...
                    dlc_pos(:,3) <= z_upper_bound & dlc_pos(:,3) >= z_lower_bound;

            % remove if ELBOW position is outside predefined cube.
            if(strcmpi(trial_data.monkey,'Han'))
                x_upper_bound = 0; % remove if larger than this
                x_lower_bound = -30; %remove if smaller than this
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
            dlc_pos = dlc_pos - median(dlc_pos,1,'omitnan');

            keep_mask_pos = keep_mask_pos & dlc_pos(:,1) <= x_upper_bound & dlc_pos(:,1) >= x_lower_bound &...
                    dlc_pos(:,2) <= y_upper_bound & dlc_pos(:,2) >= y_lower_bound & ...
                    dlc_pos(:,3) <= z_upper_bound & dlc_pos(:,3) >= z_lower_bound;   
        end
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
    
    keep_mask_all = keep_mask_spd == 1 & keep_mask_pos == 1;
end
