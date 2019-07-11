function [out] = circular_variance(in)
    % elements are in radians
    
    n = numel(in);
    if(n <=0)
        out = nan;
    else
        out = 1-sqrt(sum(cos(in)).^2 + sum(sin(in)).^2)/n;
    end

end