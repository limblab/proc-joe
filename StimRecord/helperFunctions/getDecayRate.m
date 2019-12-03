function [ output_data ] = getDecayRate( array_data, window, bin_size )
% send in a single array_data and a bin size, this then fits with a double
% exponential and outputs the decay rate from that

    array_data = rebinArrayData(array_data,bin_size);
    
    
    % get index in binCounts from window
    window_idx = [find(array_data.binEdges{1} > window(1),1,'first'),find(array_data.binEdges{1} > window(2),1,'first')];


    
    
end

