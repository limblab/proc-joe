function [ shuffled_mat ] = shuffleMatrix( mat )
% shuffles data in a NxM matrix, does not touch entries with nan's
    shuffled_mat = mat;
    
    rows = sum(sum(~isnan(mat)));
    cols = zeros(size(rows));
    data = zeros(size(rows));
    
    counter = 1;
    
    % get all non-nan data
    for i = 1:size(mat,1)
        for j = 1:size(mat,2)
            if(~isnan(mat(i,j)))
                rows(counter) = i;
                cols(counter) = j;
                data(counter) = mat(i,j);
                counter = counter + 1;
            end
        end
    end
    
    % shuffle data
    data = data(randperm(numel(data)));
    
    % put data into shuffled mat
    for i = 1:numel(data)
        shuffled_mat(rows(i),cols(i)) = data(i);
    end
    
    
end

