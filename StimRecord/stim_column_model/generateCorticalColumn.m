function [locs] = generateCorticalColumn(column_size, num)
% generates location of num neurons in a column based on the column size
% locs is a num x 3 matrix
    locs = [rand(num,1), rand(num,1), rand(num,1)].*column_size;

    % center locs so that (0,0,0) is the middle)
    locs = locs - (column_size./2);
    
    
end

