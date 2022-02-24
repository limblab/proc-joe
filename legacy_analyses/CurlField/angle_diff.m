function [ad] = angle_diff(ang1,ang2)
    % takes difference between two angles, deals with conversions
    conver_factors = [-1080,-720,-360,0,360,720,1080];
    ang1_conv = ang1 + conver_factors;
    temp_diff = ang1_conv-ang2;
    
    [~,min_idx] = min(abs(temp_diff));
    
    ad = ang1_conv(min_idx) - ang2;

end