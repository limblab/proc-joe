function [out,td] = getAdaptationMetrics(td)
    % computes time to takeoff angle relative to target for every trial and outputs
    % that in a struct
    
    takeoff_angle = zeros(numel(td),1);
    
    for tr = 1:numel(td)
        init_pos = td(tr).pos(td(tr).idx_movement_on+0,:);
        end_pos = td(tr).pos(td(tr).idx_movement_on+20,:);
        tgt_dir = td(tr).tgtDir;
        takeoff_angle(tr) = angle_diff(atan2(end_pos(2)-init_pos(2),end_pos(1)-init_pos(1))*180/pi, tgt_dir);
        
    end
    
    

    out.takeoff_angle = takeoff_angle;

end