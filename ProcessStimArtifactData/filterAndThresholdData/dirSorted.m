function [ fileListSorted ] = dirSorted( inputString )
% this function is a wrapper to dir, but sorts the files in numerical
% order. Dir does 1,10,11,...,2,21,22,...,3,.... This function does
% 1,2,3,4...,10,11,12,...,20,21,... because that lines up with the time the
% files were taken

fileList = dir(inputString); % get fileList.
fileNumbers = zeros(size(fileList,1),1);
if(size(fileList,1) > 1)
    % get all file numbers. These are after an underscore
    for fnum = 1:size(fileList,1)
        [~,fname,~] = fileparts(fileList(fnum).name);
        underscoreIdx = strfind(fname,'_');
        for uIdx = 1:numel(underscoreIdx)
            if(uIdx == numel(underscoreIdx))
                [num,status] = str2num(fname(underscoreIdx(uIdx)+1:end));
            else
                [num,status] = str2num(fname(underscoreIdx(uIdx)+1:underscoreIdx(uIdx+1)-1));
            end
            
            if(status ~= 0 && num < 1000) % check for dates
                fileNumbers(fnum) = num;
            end
        end
    end

    % sort fileNumbers, apply sort to fileList, return sorted list
    [~,sortIdx] = sort(fileNumbers);
    fileListSorted = fileList(sortIdx);
else
    fileListSorted = fileList;
end


end

