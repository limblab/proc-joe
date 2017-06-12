function [meanCounts, plusMinus] = bootstrapConfidenceInterval(cds, nn, sequenceTimes, eventTimes, binSize)
% bootstraps for 95% conf interval around mean (also provides mean) for the 
% spike counts in a 50ms bin 

cdsSpikes = cds.units(nn).spikes.ts; 
binEdges = 0:binSize:max(cds.units(nn).spikes.ts);
[binCounts, binEdges] = histcounts(cdsSpikes,binEdges);
binsRemove = 3;

% remove bins around stimulations
for et = 2:numel(eventTimes)
    % find bin closest to eventTimes(et), remove 6 bins around it or
    % something
    binIdx = floor(eventTimes(et)/binSize)-7*(et-1);
    binCounts = [binCounts(1:min(numel(binCounts)-2,binIdx-binsRemove)),binCounts(min(binIdx+binsRemove,numel(binCounts)-1):end)];
end

% remove bins before stimulation
binsPreStim = ceil(eventTimes(1)/binSize) - 8;
binCounts = binCounts(binsPreStim:end);
numSims = 10000;

% sample
% binCounts = repmat(binCounts,numSims,1);
samples = binCounts(ceil(rand(numSims,size(binCounts,2)).*size(binCounts,2))); % sample from spikes    
binMeans = mean(samples,2);

binMeans = sort(binMeans);
meanCounts = 0;
plusMinus = [0,0];
plusMinus(2) = binMeans(ceil(numSims*0.99));


% binCountsAll = sort(binCountsAll,1);
% binCountsMax = binCountsAll(ceil(numSims*0.95),:);
% plusMinus(2) = max(binCountsMax);
% % stimStartTime = eventTimes(1);
% spikes = cds.units(nn).spikes.ts; % data to sample from
% % binEdges = [binSize, 2*binSize];
% % binEdges = [-1000,binEdges,1000];
% % binEdges = 0:binSize:stimStartTime;
% % binEdges = stimStartTime:binSize:cds.meta.duration;
% sims = 1000;
% meanBinCount = zeros(sims,1);
% 
% meanCounts = 0;
% rng('shuffle');
% for i = 1:sims
%     samples = spikes(ceil(rand(size(spikes))*size(spikes,1))); 
%     sp = [];
%     spikeMask = [samples < eventTimes(1)];
%     sp = samples(spikeMask);
% %     for j = 1:length(eventTimes)
% %         if(j==length(eventTimes))
% %             spikeMask = (samples  >= eventTimes(j));
% %             sp = [sp; samples(spikeMask) - eventTimes(j)];
% %         else
% %             spikeMask = (samples  >= eventTimes(j) & ...
% %                 samples < eventTimes(j+1));
% %             sp = [sp; samples(spikeMask) - eventTimes(j)];
% %         end
% %     end
%     [bC,] = histcounts(sp,binEdges);
% %     time = cds.meta.duration - eventTimes(1) - numel(eventTimes)*binSize;
% %     meanBinCount(i) = (bC(1) + bC(3))/time; % not the bin we care about 
%     meanBinCount(i) = mean(bC);
% end
% 
% meanBinCount = sort(meanBinCount);
% plusMinus = [meanBinCount(floor(sims*0.025)), meanBinCount(ceil(sims*0.975))];
% % meanCounts = meanBinCount(floor(sims*0.5));
% % plusMinus = 2*meanCounts - plusMinus;
end

