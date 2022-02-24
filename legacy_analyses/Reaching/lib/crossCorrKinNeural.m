function [output_data] = crossCorrKinNeural(td_list)

    % compute cross correlation (with lags) between neural data and arm
    % kinematics
    % store correlation at the current lag
    % store best lag
    % store correlation at best lag
    
    r_best = zeros(size(td_list{1}.LeftS1_FR,2),numel(td_list));
    r_curr = zeros(size(r_best));
    lag_best = zeros(size(r_best));
    
    max_lag = 25; 
    
    for i_td = 1:numel(td_list)      
        
        markername = 'elbow1';
        dlc_idx = [find((strcmpi(td_list{i_td}.dlc_pos_names,[markername,'_x']))),find((strcmpi(td_list{1}.dlc_pos_names,[markername,'_y'])))];

        marker_vel = td_list{i_td}.dlc_vel(:,dlc_idx);
        marker_speed = sqrt(sum(marker_vel.^2,2));
        
        FR = td_list{i_td}.LeftS1_FR;
        
        for i_unit = 1:size(FR,2)
            [r_temp,lags_temp] = xcorr(FR(:,i_unit), marker_speed, max_lag,'coeff');
            
            % store vals
            r_curr(i_unit,i_td) = r_temp(max_lag+1);
            [r_best(i_unit,i_td), max_idx] = max(r_temp);
            lag_best(i_unit,i_td) = lags_temp(max_idx);
            
        end
    end

    % output vars
    output_data.r_curr = r_curr;
    output_data.r_best = r_best;
    output_data.lag_best = lag_best;


end