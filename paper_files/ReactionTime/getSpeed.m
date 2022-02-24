function [td] = getSpeed(td)

    for tr = 1:numel(td)
        td(tr).speed = sqrt(sum(td(tr).vel.^2,2));
        
    end


end