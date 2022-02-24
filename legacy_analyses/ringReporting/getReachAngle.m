function [ out ] = getReachAngle( td,trial_idx,center_pos )
% center_pos = [x,y] of center target. td and trial_idx within td

    out = atan2(td(trial_idx).pos(td(trial_idx).idx_otHoldTime,2)-center_pos(2),...
            td(trial_idx).pos(td(trial_idx).idx_otHoldTime,1)-center_pos(1));


end

