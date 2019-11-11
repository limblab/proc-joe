function [b_out] = getSlopeDuringTrain(array_data,condition_compute,slope_window,bin_size)

    if(isempty(condition_compute))
        condition_compute = 1:1:numel(array_data.binCounts);
    end
    
    b_out = zeros(2,numel(condition_compute));
    condition_counter = 1;
    for condition = condition_compute
        array_data = rebinArrayData({array_data},bin_size);
        array_data = array_data{1};


        slope_window_idx = [find(array_data.binEdges{condition} > slope_window(1),1,'first'),find(array_data.binEdges{condition} > slope_window(2),1,'first')];
        slope_y_vals = array_data.binCounts{condition}(slope_window_idx(1):slope_window_idx(2))/array_data.num_stims(condition)/(bin_size)*1000;
        slope_x_vals = array_data.binEdges{condition}(slope_window_idx(1):slope_window_idx(2));
        b_out(:,condition_counter) = [ones(numel(slope_y_vals),1), slope_x_vals']\slope_y_vals';
        condition_counter = condition_counter + 1;
    end

end