function [meanCounts, plusMinus] = bootstrapConfidenceInterval(cds, nn, stimStartTime,binSize)
% bootstraps for 95% conf interval around mean (also provides mean) for the 
% spike counts in a 50ms bin 

spikes = cds.units(nn).spikes.ts; % data to sample from
binEdges = 0:0.05:stimStartTime;
sims = 1000;
meanBinCount = zeros(sims,1);
for i = 1:sims
    spikeSamples = spikes(ceil(rand(size(spikes))*size(spikes,1))); 
    spikeMask = [cds.units(nn).spikes.ts < stimStartTime]; 
    spikeSamples = spikeSamples(spikeMask);
    [bC,] = histcounts(spikeSamples,binEdges);
    meanBinCount(i) = mean(bC); 
end

meanBinCount = sort(meanBinCount);
plusMinus = [meanBinCount(floor(sims*0.025)), meanBinCount(ceil(sims*0.975))];
meanCounts = meanBinCount(floor(sims*0.5));
% plusMinus = 2*meanCounts - plusMinus;
end

