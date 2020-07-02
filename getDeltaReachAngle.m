function [output_data] = getDeltaReachAngle(td,params)

    output_data = [];

    mean_go_cue_stim_time = ceil(mean([td.idx_stimTime] - [td.idx_goCueTime],'omitnan'));
    
    [reach_angles, tgt_dir_list, mean_angle, std_angle] = getReachAngleDistribution(td,mean_go_cue_stim_time);
    tgt_dirs = unique(tgt_dir_list);
    
    td_stim = td(~isnan([td.idx_stimTime]));
    
    delta_angle = nan(numel(td_stim),1);
    stim_code = nan(numel(td_stim),1);
    tgt_dir = nan(numel(td_stim),1);
    correct_latency = nan(numel(td_stim),1);
    
    for i_trial = 1:numel(td_stim)
        stim_code(i_trial) = td_stim(i_trial).stimCode;
        tgt_dir(i_trial) = td_stim(i_trial).target_direction;
        tgt_dir_idx = find(tgt_dirs == td_stim(i_trial).target_direction);
        
        % get stim idx
        if(~isnan(stim_code(i_trial)))
            idx_stim = td_stim(i_trial).idx_stimTime;
        else
            idx_stim =  td_stim(i_trial).idx_goCueTime + mean_go_cue_stim_time;
        end
        
        % get correction latency relative to idx_stim
        % find when reach_angle is greater than mean_angle*std_angle
        data = td_stim(i_trial).vel(idx_stim-1:idx_stim+ceil(0.3/td_stim(i_trial).bin_size),:);
        trial_angles = atan2(data(:,2),data(:,1));
        trial_angles = trial_angles(2:end) - trial_angles(1);
        
        idx_traj_change = find(circ_dist(trial_angles,mean_angle(tgt_dir_idx)) > 2*std_angle(tgt_dir_idx));
        
        if(~isempty(idx_traj_change) & idx_traj_change(1) > 0.12/td_stim(i_trial).bin_size)
            % make sure idx_traj_change does not occur at edge of
            % workspace
            
            
            delta_angle(i_trial) = circ_mean(circ_dist(trial_angles(idx_traj_change),mean_angle(tgt_dir_idx)));
            correct_latency(i_trial) = idx_traj_change(1)*td_stim(i_trial).bin_size;
        end
        
    end
    
    % format output_data
    
    output_data.delta_angle = delta_angle;
    output_data.stim_code = stim_code;
    output_data.tgt_dir = tgt_dir;
    output_data.correct_latency = correct_latency;
    
end


function [reach_angles, tgt_dir_list, mean_angle, std_angle] = getReachAngleDistribution(td,go_cue_offset)

    reach_angles = []; tgt_dir_list = []; mean_angle = []; std_angle = [];

    td_no_stim = td(isnan([td.idx_stimTime]) | [td.idx_stimTime] < [td.idx_goCueTime]);
    
    max_offset = go_cue_offset + ceil(0.2/td_no_stim(1).bin_size);
    
    for i_trial = 1:numel(td_no_stim)
        data = td_no_stim(i_trial).vel(td_no_stim(i_trial).idx_goCueTime+go_cue_offset-1:td_no_stim(i_trial).idx_goCueTime+max_offset,:);
    
        
        trial_angles = atan2(data(:,2),data(:,1));
        trial_angles = trial_angles(2:end)-trial_angles(1);
        
        reach_angles(end+1:end+numel(trial_angles)) = trial_angles;
        tgt_dir_list(end+1:end+numel(trial_angles)) = td_no_stim(i_trial).target_direction;
    end
    
    tgt_dirs = unique(tgt_dir_list);
    
    for i_tgt = 1:numel(tgt_dirs)
        mean_angle(end+1) = circ_mean(reach_angles(tgt_dir_list==tgt_dirs(i_tgt))');
        % circ_var provides 1-R, we want sqrt(ln(-2*R))
        std_angle(end+1) = sqrt(-2*log((1-circ_var(reach_angles(tgt_dir_list==tgt_dirs(i_tgt))')))); 
    end

end
