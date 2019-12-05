function [ output_data ] = getDecayRate( array_data, window, bin_size )
% send in a single array_data and a bin size, this then fits with a double
% exponential and outputs the decay rate from that

    array_data = rebinArrayData(array_data,bin_size);
    
    
    % get index in binCounts from window
    window_idx = [find(array_data.binEdges{1} > window(1),1,'first'),find(array_data.binEdges{1} > window(2),1,'first')];

    % for each condition, fit data
    fits = cell(size(array_data.binCounts));
    gof = cell(size(fits));
    param_list = zeros(numel(fits),2);
    
    for condition = 1:numel(fits)
        x_data = array_data.binEdges{condition}(window_idx(1):window_idx(2))'/1000; % x_data in seconds
        y_data = array_data.binCounts{condition}(window_idx(1):window_idx(2))';
        
        % fit data
        [fits{condition},gof{condition}] = fit(x_data,y_data,'a*exp(-b*x)','StartPoint',[10,1]);
        
        param_list(condition,:) = [fits{condition}.a, fits{condition}.b];
        
    end
    
    
    
    % format outputs
    output_data.fits = fits;
    output_data.gof = gof;
    output_data.param_list = param_list;
    
end

