function [out] = circular_mean(in)
    % elements are in radians
    
    n = numel(in);
    if(n <=0)
        out = nan;
    else
        out = atan2(1/n*sum(sin(in)),1/n*sum(cos(in)));
    end

end