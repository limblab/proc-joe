function [prop_data] = getModelCellTypeProportions()

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