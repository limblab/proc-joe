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
        markername = 'hand2';
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
        %dlc_rot(:,1:2)
        %handle_pos
        %sqrt(sum((dlc_rot(:,1:2)-handle_pos).^2,2))
        keep_mask = sqrt(sum((dlc_rot(:,1:2)-handle_pos).^2,2)) < 1000; % 3.5 was chosen again by looking at a histogram of distances for an example dataset
        
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
        
        if(strcmpi(trial_data.monkey,'Han'))
            % we have already moved bad tracking, just remove if hand is on
            % the face plate. We have predefined (x,y,z) bounds for this
            
            x_bound = 16; % remove if larger than this, and others
            x_lower_bound = -30; %remove if smaller than this
            y_bound = 3; % remove if smaller than this, and others
            z_bound = 4.5; % remove if larger than this, and others
        else          
            x_bound = 20; % remove if larger than this, and others
            x_lower_bound = -30; %remove if smaller than this
            y_bound = 10; % remove if smaller than this, and others
            z_bound = 5; % remove if larger than this, and others
            
            markername = 'hand2';
            dlc_idx = [find(strcmpi(trial_data.dlc_pos_names,[markername,'_x'])),find(strcmpi(trial_data.dlc_pos_names,[markername,'_y'])),...
                find(strcmpi(trial_data.dlc_pos_names,[markername,'_z']))];
            
            keep_mask = ~(trial_data.dlc_pos(:,dlc_idx(1)) > x_bound & ...
                trial_data.dlc_pos(:,dlc_idx(2)) < y_bound & ...
                trial_data.dlc_pos(:,dlc_idx(3)) > z_bound);
            %disp(keep_mask == 0)
            marker_spd = calcSpeed(trial_data.dlc_pos(:,dlc_idx(1)), trial_data.dlc_pos(:,dlc_idx(2)), trial_data.dlc_pos(:,dlc_idx(3)),trial_data.bin_size);
            %plot(marker_spd)
            %disp(mean(marker_spd))
            
            %The second pointer of this two-pointer function. When the
            %average speed of N (determined by judgingLength parameter
            %below) bins is larger than a threshold (determined by 
            %inExpPhaseThreshold), set this switch to 1, else, set it to 0.
            %If the expPhaseSwitch is switched from 1 to 0, record these
            %bins as "1" in the keep_mask array. Else if the
            %expPhaseSwitch is switched from 0 to 1, record these bins
            %(all the consecutive "0" bins) as "0" in the keep_mask array.
            expPhaseSwitch = 0;
            
            %This parameter determines how many bins of speed data
            %the code is going to take average of. This averaged speed data
            %will be used to judge whether this bin is in experiment
            %phase or not.
            judgingLength = 10; %bins, currently each bin is 0.05s.
            
            %If for each bin, the average value of its upcoming N frames
            %(determined by judgingLength) is higher than this threshold,
            %consider it in the experiment phase already.
            inExpPhaseThreshold = mean(marker_spd);
            
            %disp(size(marker_spd))
            %disp(marker_spd(1:5))
            %disp(marker_spd(1:5,1))
            
            binNumPointer = 1;
            
            %Preliminarily predict if any data POINTS are in experiment
            %phase
            for i = 1:length(marker_spd)-judgingLength
                curr_mean_spd = mean(marker_spd(i:i+judgingLength));
                if curr_mean_spd > inExpPhaseThreshold && expPhaseSwitch == 0
                    expPhaseSwitch = 1;
                    keep_mask(binNumPointer:i) = 0;
                    binNumPointer = i;
                end
                if curr_mean_spd < inExpPhaseThreshold && expPhaseSwitch == 1
                    expPhaseSwitch = 0;
                    keep_mask(binNumPointer:i) = 1;
                    binNumPointer = i;
                end
            end
            
            %I'm not sure how to properly explain this variable...
            %Basically what I want to do is, skim through the datset, and
            %if a point I and its continuous N points (which is this
            %variable) has an average value over M (another variable), we
            %say that ALL these points (I:I+N) are considered IN EXPERIMENT
            %PHASE
            scanLen = 40; %frames/bins
            expPhaseThreshold = 0.15; %If 60% of the points in this section were considered as 1
            
            i = 1;
            %COMBINE these points into large experiment-phase segments
            %for i = 1:length(keep_mask)
            while i < length(keep_mask) - scanLen
                if mean(keep_mask(i:i+scanLen)) > expPhaseThreshold
                    keep_mask(i:i+scanLen) = 1;
                    i = i + scanLen;
                else
                    keep_mask(i:i+scanLen) = 0;
                    i = i + scanLen;
                end
                i = i + 1;
            end

            
            %Maximum speed threshold to remove data points where the
            %markers are clearly way over speed
            maxSpeedThres = 2
            
            keep_mask(marker_spd > maxSpeedThres) = 0;
            
            
            
            
            
            
            %Plot the actual speed data and the predicted result to see if
            %they have aligned well
            X_axis = 0:trial_data.bin_size:(length(marker_spd))*trial_data.bin_size;
            length(X_axis)
            length(marker_spd)
            figure;
            plot(X_axis(1:length(marker_spd)),marker_spd);
            hold on
            k = plot(X_axis(1:length(keep_mask)),keep_mask,'LineWidth',1);
            hold off 
            %alpha(k,.2);
            %xlim([0, length(marker_spd)*trial_data.bin_size]);
            legend('Actual Marker Speed','Predicted Mask');
            
            
            
            
            
            
            

            % remove points
            td_fields = fieldnames(trial_data);
            pos_len = length(trial_data.dlc_pos);
            for i_field = 1:numel(td_fields)
                if(length(trial_data.(td_fields{i_field})) == pos_len)
                    trial_data.(td_fields{i_field}) = trial_data.(td_fields{i_field})(keep_mask==1,:);
                end
            end
            trial_data.trial_idx = find(keep_mask);
            
        %else
        %    warning('getExperimentPhase is not implemented for this task and monkey. Nothing was done');
        end
        
    end %end of either process the RT2D task or the RT3D task to experiment-only data
    
    function [speed] = calcSpeed(pos_x, pos_y, pos_z, bin_size)
        vel_x = pos_x(2:length(pos_x)) - pos_x(1:length(pos_x)-1);
        vel_x = [0.0; vel_x];
        vel_y = pos_y(2:length(pos_y)) - pos_y(1:length(pos_y)-1);
        vel_y = [0.0; vel_y];
        vel_z = pos_z(2:length(pos_z)) - pos_z(1:length(pos_z)-1);
        vel_z = [0.0; vel_z];
        speed = sqrt(vel_x.^2 + vel_y.^2 + vel_z.^2) / 100 / bin_size;
    end
    




end