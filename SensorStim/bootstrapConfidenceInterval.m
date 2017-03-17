function [meanCounts, plusMinus] = bootstrapConfidenceInterval(cds, nn, stimStartTime,binSize)
% bootstraps for 95% conf interval around mean (also provides mean) for the 
% spike counts in a 50ms bin 

spikes = cds.units(nn).spikes.ts; % data to sample from
binEdges = 0:binSize:stimStartTime;
% binEdges = stimStartTime:binSize:cds.meta.duration;
sims = 1000;
spikeMask = [cds.units(nn).spikes.ts < stimStartTime]; 
sp = spikes(spikeMask);
[binCounts, binEdges] = histcounts(sp,binEdges);
meanBinCount = zeros(sims,1);

meanCounts = mean(binCounts);
rng('shuffle');
for i = 1:sims
    samples = spikes(ceil(rand(size(spikes))*size(spikes,1))); 
    spikeMask = [samples < stimStartTime]; 
    samples = samples(spikeMask);
    [bC,] = histcounts(samples,binEdges);
    meanBinCount(i) = mean(bC); 
end

meanBinCount = sort(meanBinCount);
plusMinus = [meanBinCount(floor(sims*0.025)), meanBinCount(ceil(sims*0.975))];
% meanCounts = meanBinCount(floor(sims*0.5));
% plusMinus = 2*meanCounts - plusMinus;
end

