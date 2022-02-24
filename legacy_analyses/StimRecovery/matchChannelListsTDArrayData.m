function [td_move, arrayData] = matchChannelListsTDArrayData(td_move, arrayData,array_name)

    td_chan_list = td_move(1).([array_name,'_unit_guide'])(:,1);
    td_id_list = td_move(1).([array_name,'_unit_guide'])(:,2);
    
    arr_chan_list = [];
    arr_id_list = [];
    td_keep_mask = ones(numel(td_chan_list),1);
    arr_keep_mask = ones(numel(arrayData),1);
    for arr_idx = 1:numel(arrayData)
        arr_chan_list(arr_idx) = arrayData{arr_idx}.CHAN_REC;
        arr_id_list(arr_idx) = arrayData{arr_idx}.ID;
        
    end

    % need those two lists to match....so either remove from one or keep in
    % both
    common_chan_list = union(arr_chan_list,td_chan_list);
    for i_chan = 1:numel(common_chan_list)
        % go through each ID for this channel and see if they exist in both
        % arrays
        common_ids = union(arr_id_list(arr_chan_list == common_chan_list(i_chan)),...
            td_id_list(td_chan_list == common_chan_list(i_chan)));
        for i_id = 1:numel(common_ids)
            % find idx for both
            td_idx = find(td_chan_list == common_chan_list(i_chan) & td_id_list == common_ids(i_id));
            arr_idx = find(arr_chan_list == common_chan_list(i_chan) & arr_id_list == common_ids(i_id));
            if(isempty(td_idx))
                arr_keep_mask(arr_idx) = 0;
            end
            if(isempty(arr_idx))
                td_keep_mask(td_idx) = 0;
            end
        end
    end
    
    % adjust trial data and arrayData
    for i_trial = 1:numel(td_move)
        chan_mask = any(td_move(i_trial).([array_name,'_unit_guide'])(:,1) + 1000*td_move(i_trial).([array_name,'_unit_guide'])(:,2) == arr_chan_list,2);
        
        td_move(i_trial).([array_name,'_spikes']) = td_move(i_trial).([array_name,'_spikes'])(:,td_keep_mask==1);
        td_move(i_trial).([array_name,'_unit_guide']) = td_move(i_trial).([array_name,'_unit_guide'])(td_keep_mask==1,:);
    end
    arrayData = arrayData(arr_keep_mask==1);



end