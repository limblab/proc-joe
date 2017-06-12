function [meanCounts, plusMinus] = bootstrapConfidenceIntervalSpindle(cds, nn, sequenceTimes, eventTimes, binSize, stimState)
% bootstraps for 95% conf interval around mean (also provides mean) for the 
% spike counts in a 50ms bin 

timeDiff = 1000;
stimStart = eventTimes(1);
cdsSpikes = cds.units(nn).spikes.ts(cds.units(nn).spikes.ts(:) < stimStart); % data to sample from
numSims = 1000;
% initialize bins and set edges
binEdges = 0:binSize:stimStart;
[binCounts, binEdges] = histcounts(cdsSpikes,binEdges);

meanSims = zeros(numSims,1);
for sim = 1:numSims
    samples = binCounts(ceil(rand(size(binCounts))*size(binCounts,2))); % sample from spikes
    meanSims(sim,1) = mean(samples);
end
meanSims = sort(meanSims);
meanCounts = 0;
plusMinus = [0,0];
plusMinus(2) = meanSims(ceil(numSims)*0.99);
if(plusMinus(2) == 0)
    plusMinus(2) = 1000;
end
% cdsSpikes = cds.units(nn).spikes.ts; % data to sample from
% preTime = max(eventTimes(:) - sequenceTimes(:,1));
% postTime = max(sequenceTimes(:,2)-eventTimes(:));
% numBinsPre = ceil(preTime/binSize);
% numBinsPost = ceil(postTime/binSize);
% numSims = 1000;
% % initialize bins and set edges
% binEdges = [0];
% for i = 1:numBinsPre
%     binEdges = [-i*binSize binEdges];
% end
% for i = 1:numBinsPost
%     binEdges = [binEdges i*binSize];
% end
% 
% stimSampRate = 1;
% if(~isempty(cds.analog))
%     stimSampRate = 1/(cds.analog{1,1}.t(2) - cds.analog{1,1}.t(1));
% else
%     stimSampRate = 1/(cds.lfp.t(2) - cds.lfp.t(1));
% end
% aveStimTime = (sum(stimState == 1))/numel(eventTimes)/stimSampRate;
% 
% binCare = numBinsPre + [1:ceil(aveStimTime/binSize) + 2];
% binCountsAll = zeros(numSims, numBinsPre+numBinsPost-numel(binCare));
% 
% for sim = 1:numSims
%     samples = cdsSpikes(ceil(rand(size(cdsSpikes))*size(cdsSpikes,1))); % sample from spikes
%     
%     % adjust timing of spikes based on eventTimes
%     spikes = [];
%     for i = 1:length(eventTimes)
%         spikeMask = (samples  > eventTimes(i) - preTime & ...
%             samples < eventTimes(i) + postTime);
%         spikes = [spikes; samples(spikeMask) - eventTimes(i)];
%     end
% 
%     [binCounts, binEdges] = histcounts(spikes,binEdges);
%     binCounts = binCounts/numel(eventTimes);
%     binCountsAll(sim,:) = [binCounts(1:binCare(1)-1),binCounts(binCare(end)+1:end)];
% end
% 
% meanCounts = 0;
% plusMinus = [0,0];
% % binCountsMax = max(binCountsAll,[],2);
% % binCountsMax = sort(binCountsMax);
% % plusMinus(2) = binCountsMax(ceil(numSims*0.95));
% binCountsAll = sort(binCountsAll,1);
% binCountsMax = binCountsAll(ceil(numSims*0.95),:);
% plusMinus(2) = max(binCountsMax);


end

