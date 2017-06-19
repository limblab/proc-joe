function [ out ] = getRawDataIdx(unitTs,unitChan,rawTs,rawChan)
% unitTs and unitChan are single numbers, rawTs and rawChan are all of them
out = -1;
rawIdx = find(unitTs == rawTs);

if(numel(rawIdx) > 1)
    for i = 1:numel(rawIdx)
        if(rawChan(rawIdx(i)) == unitChan)
            out = rawIdx(i);
        end
    end
elseif(numel(rawIdx) ~= 0)
    out = rawIdx;
end


end

