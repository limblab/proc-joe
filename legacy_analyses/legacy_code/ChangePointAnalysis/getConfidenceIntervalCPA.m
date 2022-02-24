function [ confInter ] = getConfidenceIntervalCPA(data,cutoff)
confInter = [0,0];
% find all intersections
[~,maxIdx] = max(data);

dataLow = data(1:maxIdx-1);
dataHigh = data(maxIdx+1:end);

dataLow(dataLow > cutoff) = 100*cutoff;
dataHigh(dataHigh > cutoff) = 100*cutoff;

[~,confInter(1)] = min(abs(dataLow-cutoff));
[~,temp] = min(abs(dataHigh-cutoff));
confInter(2) = temp + numel(dataLow) + 1;

end

