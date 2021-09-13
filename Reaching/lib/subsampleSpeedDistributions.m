function [td_sub] = subsampleSpeedDistributions(td_list,task_list)

    
    % for every data point in the 3D task, try and find a data point in the
    % 2D task with a similar speed. If within some bound (5 cm/s?)
    
    speed_bound = 10;
    
    td_idx_3d = find(strcmpi(task_list,'RT3D'));
    td_idx_2d = find(strcmpi(task_list,'RT'),1,'first');
    
    td_3d = td_list{td_idx_3d};
    td_2d = td_list{td_idx_2d};
    
    % get hand speed for each trial data
    markername = 'hand2';
    dlc_idx = [find(strcmpi(td_3d.dlc_pos_names,[markername,'_x'])),find(strcmpi(td_3d.dlc_pos_names,[markername,'_y'])),...
        find(strcmpi(td_3d.dlc_pos_names,[markername,'_z']))];

    % go through data set and find long periods of low speed, mark as
    % non-experiment
    spd_3d = sqrt(sum(td_3d.dlc_vel(:,dlc_idx).^2,2));
    spd_2d = sqrt(sum(td_2d.dlc_vel(:,dlc_idx).^2,2));
    
    keep_2d = zeros(size(spd_2d));
    keep_3d = zeros(size(spd_3d));
    
    order_3d = randperm(numel(keep_3d));
    
    for i_3d = 1:numel(order_3d)
        idx_3d = order_3d(i_3d);
        % search for speed within speed bound
        valid_point = abs(spd_3d(idx_3d) - spd_2d) < speed_bound & keep_2d == 0;
        
        if(sum(valid_point) > 0) % randomly choose one
            keep_3d(idx_3d) = 1;
            
            poss_2d_idx = find(valid_point);
            valid_2d_idx = datasample(poss_2d_idx,1);
            
            keep_2d(valid_2d_idx) = 1;
        
        %else % no valid points, keep_3d remains as 0
        end
    end
    
    
    % keep wanted data
    td_names = fieldnames(td_3d);
    field_size = size(td_3d.dlc_pos,1);
    for i_name = 1:numel(td_names)
        if(field_size == size(td_3d.(td_names{i_name}),1))
            td_3d.(td_names{i_name})(~keep_3d,:) = [];
        end
    end
    
    td_names = fieldnames(td_2d);
    field_size = size(td_2d.dlc_pos,1);
    for i_name = 1:numel(td_names)
        if(field_size == size(td_2d.(td_names{i_name}),1))
            td_2d.(td_names{i_name})(~keep_2d,:) = [];
        end
    end
    
    % package outputs
    td_sub = {};
    td_sub{td_idx_3d} = td_3d;
    td_sub{td_idx_2d} = td_2d;





end