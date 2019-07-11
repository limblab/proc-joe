function [out] = mean_cosine_error(data,pred)
    % elements are in radians, if 2d array is put in, then acts upon second
    % dimension (returns an nx1 array if an nxm array is sent in)
    if(any(size(data) == 1) && iscolumn(data))
        data = data';
    end
    
    n = numel(data);
    if(n <=0)
        out = nan;
    else
        out = 1-1/(size(data,2))*sum(cos(data-pred),2);
    end

end