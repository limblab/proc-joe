function [td] = setOriginAsShoulder(td,is_fixed_origin)
    

    dlc_idx = [find((strcmpi(td.dlc_pos_names,'shoulder_x'))),...
        find((strcmpi(td.dlc_pos_names,'shoulder_y'))),...
        find((strcmpi(td.dlc_pos_names,'shoulder_z')))];

    if(is_fixed_origin)
        
        dlc_pos = td.dlc_pos(:,dlc_idx);
        dlc_pos = dlc_pos(~any(isnan(dlc_pos),2),:);
        
        origin = dlc_pos(1,:);
        
    else
        dlc_pos = td.dlc_pos(:,dlc_idx);
        isnan_mask = any(isnan(dlc_pos),2);
        dlc_pos = dlc_pos(~isnan_mask,:);
        
        origin_nan = dlc_pos(1,:);
        td.dlc_pos(isnan_mask==1,dlc_idx(1)) = origin_nan(1);
        td.dlc_pos(isnan_mask==1,dlc_idx(2)) = origin_nan(2);
        td.dlc_pos(isnan_mask==1,dlc_idx(3)) = origin_nan(3);
        
        origin = td.dlc_pos(:,dlc_idx);
    end

    x_idx = find(~cellfun(@isempty,strfind(td.dlc_pos_names,'_x')));
    y_idx = find(~cellfun(@isempty,strfind(td.dlc_pos_names,'_y')));
    z_idx = find(~cellfun(@isempty,strfind(td.dlc_pos_names,'_z')));

    td.dlc_pos(:,x_idx) = td.dlc_pos(:,x_idx)-origin(:,1);
    td.dlc_pos(:,y_idx) = td.dlc_pos(:,y_idx)-origin(:,2);
    td.dlc_pos(:,z_idx) = td.dlc_pos(:,z_idx)-origin(:,3);
   


end