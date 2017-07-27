function [ out ] = getRawDataIdx(unitTs,unitChan,rawTs,rawChan)
% unitTs and unitChan are single numbers, rawTs and rawChan are all of them
out = -1;
unitTs = round(unitTs/0.001)*0.001;
rawTs = round(rawTs/0.001)*0.001;

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

