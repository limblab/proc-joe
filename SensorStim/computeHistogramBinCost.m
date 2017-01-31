function [ cost ] = computeHistogramBinCost( spikes, binSize, dt )
%% computes the cost (C = (2*k-v)/(bs^2) given the binCounts and binSize
numBins = floor(dt/binSize);
[binCounts, binEdges] = histcounts(spikes,numBins);
meanBins = mean(binCounts);
varBins = sum((binCounts-meanBins).^2)/numel(binCounts);
cost = (2*meanBins - varBins)/(binSize.^2);

end

