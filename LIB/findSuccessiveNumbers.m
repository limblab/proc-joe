function [ idx_list ] = findSuccessiveNumbers( vec, num_needed )
% finds index of successive numbers in a vector,need at least num in a row 
% assumes vec is a list of numbers (n x 1 or 1 x n)

% idx_list = [idx_start(1), idx_end(1); idx_start(2), idx_end(2)...] for
% every group of successive numbers

    idx_list = [];
    is_suc = 0;
    idx_start = [];
    num_in_suc = 1; % starts as 1 bc current item counts
    
    for v = 1:(numel(vec)-1)
        if(vec(v) == vec(v+1) - 1)
            is_suc = 1;
            num_in_suc = num_in_suc + 1;
            if(isempty(idx_start))
                idx_start = v;
            end
        elseif(is_suc) % successive list
            if(num_in_suc >= num_needed) % store idx_start and v
                idx_list(end+1,:) = [idx_start,v];                 
            end
            % start counter over
            num_in_suc = 1; is_suc = 0; idx_start = [];
        end
    end



end

