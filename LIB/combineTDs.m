function [td_comb] = combineTDs(td_one,td_two)

    % combines trial datas with different neural fields (lengths) -- aligns
    % based on electrode ID
    td_comb = [];
    if(isempty(td_one))
        td_comb = td_two;
        return; % nothing to combine
    end
    
    if(isempty(td_two))
        td_comb = td_one;
        return; % nothing to combine
    end
    
    % get field name -- find _unit_guide
    fields = fieldnames(td_one);
    unit_guide_idx = -1;
    for i_field = 1:numel(fields)
        if(~isempty(strfind(fields{i_field},'unit_guide')))
            unit_guide_idx = i_field;
        end
    end
    array_name = fields{unit_guide_idx}(1:end-11);
    
    % get electrode and unit data for both trial datas
    neural_data_one = td_one.([array_name,'_unit_guide']);
    neural_data_two = td_two.([array_name,'_unit_guide']);
    
    neural_data_union = union(neural_data_one,neural_data_two,'rows');
    
    % put each trial of td_one in td_comb
    for i_trial = 1:numel(td_one)
        td_comb = appendRowInTD(td_comb,td_one(i_trial),neural_data_union,array_name);
    end
    
    % put each trial of td_two in td_comb
    for i_trial = 1:numel(td_two)
        td_comb = appendRowInTD(td_comb,td_two(i_trial),neural_data_union,array_name);
    end
    
end



function [td_master] = appendRowInTD(td_master, td, neural_data_id, array_name)
    
    orig_unit_guide = td.([array_name,'_unit_guide']);
    % update neural fields to match neural_data_id
    spike_data = [];
    unit_guide = neural_data_id;
    ts = {};
    % remake neural fields to match neural data id
    for i_neural = 1:size(neural_data_id,1)
        % look for field with same electrode and unit as neural_data_id
        same_elec = find(neural_data_id(i_neural,1) == orig_unit_guide(:,1));
        same_ID = find(neural_data_id(i_neural,2) == orig_unit_guide(:,2));
        same_info = intersect(same_elec, same_ID);
        
        % append to spike_data if it exists, otherwise skip idx
        if(~isempty(same_info))
            spike_data(:,i_neural) = td.([array_name,'_spikes'])(:,same_info);
%             ts{i_neural,1} = td.([array_name,'_ts']){same_info};
        end
    end
    % update neural fields
    td.([array_name,'_unit_guide']) = unit_guide;
    td.([array_name,'_spikes']) = spike_data;
%     td.([array_name,'_ts']) = ts;
    
    % append to td_master
    td_master = [td_master,td];

end
