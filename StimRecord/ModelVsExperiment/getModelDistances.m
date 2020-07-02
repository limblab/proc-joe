function [dist_data] = getModelDistances(array_data)
    dist_data = zeros(numel(array_data),1);
    
    
    for i = 1:numel(array_data)
        dist_data(i) = norm(array_data{i}.loc,2);
    end


end