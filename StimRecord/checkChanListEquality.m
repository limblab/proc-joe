function [ out ] = checkChanListEquality( chanSent,chanList )
% chanSent is an nx1 cell array where each entry contains a list of
% channels stimulated. This tells you which channel(s) were stimulated on
% each trial

% chanList is a list of channels

% this function returns an nx1 array where each entry says if chanSent{i}
% is equal to chanList. All channels in chanList must be in chanSent for
% this function to return 1 in that array location

    sortedChanList = sort(chanList);
    out = ones(numel(chanSent),1); % assume equality

    for i = 1:numel(out)
        sortedChanSent = sort(chanSent{i});
        if(numel(sortedChanSent) == numel(chanList))
            for j = 1:numel(sortedChanSent)
                if(sortedChanSent(j) ~= sortedChanList)
                    out(i) = 0; % can't be same
                end
            end
        else % not the same length, can't be same list
            out(i) = 0;
        end
        
    end

end

