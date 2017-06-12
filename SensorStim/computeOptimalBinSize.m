function [ binCounts, binEdges ] = computeOptimalBinSize(spikes, initBinSize, preTime, postTime)
% computes optimal bin size given the spikes time normalized by eventTimes
% minimize C = (2*k-v)/(binSize^2); where k is mean of bins, v = var of
% bins

% Use pseudo grad descent <- I made the name up whoops
dt = postTime-preTime;
binSize = initBinSize/10;
binSizeBest = binSize;
costCurrent = computeHistogramBinCost(spikes, binSize,dt);

for i = 1:10 % ten rounds
    costLeft = computeHistogramBinCost(spikes,binSize*(1-1/(2^i)),dt);
    costRight = computeHistogramBinCost(spikes,binSize*(1+1/(2^i)),dt);
    if(costCurrent < costLeft && costCurrent < costRight)
        costCurrent = costCurrent;
        binSizeBest = binSize;
        binSize = (0.1-0.01)*rand()+0.01;
    elseif(costLeft < costRight)
        costCurrent = costLeft;
        binSize = binSize*(1-1/(2^i));
        binSizeBest = binSize;
    else
        costCurrent = costRight;
        binSize = binSize*(1+1/(2^i));
        binSizeBest = binSize;
    end
end
numBins = floor((postTime-preTime)/binSizeBest);
numBinsPre = ceil(preTime/binSizeBest);
numBinsPost = ceil(postTime/binSizeBest);
% set bin edges
binEdges = [0];
for i = 1:numBinsPre
    binEdges = [-i*binSizeBest binEdges];
end
for i = 1:numBinsPost
    binEdges = [binEdges i*binSizeBest];
end
[binCounts, binEdges] = histcounts(spikes,binEdges);

end

