function [sampled_array_data,sampled_mask_data] = resampleModelData(array_data,mask_data)

    % this function resamples (with replacement) an array data based on the
    % cell proportions from Markram 2015 (see getCellTypeProportions below)
    
    % mask data contains diameter, cell_id, and clone num  
    
    prop_data = getCellTypeProportions();
    
    sampled_array_data = {};
    sampled_mask_data = [];
    sampled_mask_data.diam = []; sampled_mask_data.cell_id = []; sampled_mask_data.clone = [];
    % sample for each diameter
    unique_diams = unique(mask_data.diam);
    for i_diam = 1:numel(unique_diams)
        idx_diam = find(mask_data.diam == unique_diams(i_diam));
        cell_id_list = mask_data.cell_id(idx_diam);
        prop_val = repmat(prop_data(:,2)',numel(idx_diam),1);
        
        prop_idx = cell_id_list == prop_data(:,1)';
        
        prop_sample = prop_val(prop_idx);
        
        % sample based on probabilities in prop_sample
        
        idx_sample = datasample(idx_diam,numel(idx_diam),'Replace',true,'Weights',prop_sample);
        
        sampled_array_data(end+1:end+numel(idx_sample)) = array_data(idx_sample);
        sampled_mask_data.diam(end+1:end+numel(idx_sample),1) = mask_data.diam(idx_sample);
        sampled_mask_data.cell_id(end+1:end+numel(idx_sample),1) = mask_data.cell_id(idx_sample);
        sampled_mask_data.clone(end+1:end+numel(idx_sample),1) = mask_data.clone(idx_sample);
    end
    
    
    
end




function [prop_data] = getCellTypeProportions()

    % cell types are estimated from Markram 2015 fig 2,3
    
    layer_num_cells = [322,3400+4099,4654,6106,12738]; % L1, L2/3, L4, L5, L6
    prop_excite_inhib = [1,0.8,0.1,0.75,0.9]; % 
    is_excite = [1,1,0,1,1];
    prop_cell = [0.2,1,0.25,0.75,0.2];
    
    % outputs a [cell_id, proportion of population which should be of that
    % cell] matrix

    prop_data = zeros(5,2);
    prop_data(:,1) = [1:5:21];
    prop_data(:,2) = layer_num_cells.*prop_excite_inhib.*prop_cell;
    prop_data(:,2) = prop_data(:,2)/sum(prop_data(:,2));
    
end