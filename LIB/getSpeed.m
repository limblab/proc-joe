function [ td ] = getSpeed( td )
% makes a speed entry in trial data

    for t = 1:numel(td)
        td(t).speed = sqrt(td(t).vel(:,1).^2 + td(t).vel(:,2).^2);
        
    end


end

