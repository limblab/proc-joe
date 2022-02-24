function [ out ] = getRawDataIdx(unitTs,unitChan,rawTs,rawChan)
% find unitTs and unitChan in rawTs and rawChan by comparing time stamps
% and channels. out holds the idx in rawTs that unitTs and unitChan
% correspond to


% initialize out
out = -1*ones(numel(unitTs),1);

% round times to nearest millisecond
unitTs = round(unitTs/0.001)*0.001; 
rawTs = round(rawTs/0.001)*0.001;

for idx = 1:numel(unitTs)
    rawIdx = find(unitTs(idx) == rawTs);

    if(numel(rawIdx) > 1)
        for i = 1:numel(rawIdx)
            if(rawChan(rawIdx(i)) == unitChan(idx))
                out(idx) = rawIdx(i);
            end
        end
    elseif(numel(rawIdx) ~= 0)
        out(idx) = rawIdx;
    end
end


end

