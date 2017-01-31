function [ histogramOut ] = generatePESTH(cds,neuronNumber, varargin)
% possible varargin: 'useRate', 1 <- plots rate 
%                    'optimalBinSize', 1 <-calculates and plots PSTH with optimal bin size as well, 
%                    'binSize', # <- uses this binSize instead of 0.05.

% initialize variables
userBinSize = 0.05;
optimalBinSize = 0;
useRate = 0;
% deal with varargin
for i = 1:2:size(varargin,2)
    switch varargin{i}
        case 'useRate'
            useRate = varargin{i+1};
        case 'optimalBinSize'
            optimalBinSize = varargin{i+1};
        case 'binSize'
            userBinSize = varargin{i+1};
    end
end


warning('off')
funcFolder = pwd;

% generates a post stimulus time histogram from the given unit in the cds
GTOstim = 0;
if(size(cds.analog,1) == 0)
    GTOstim = 1;
end

[stimState,] = determineStimTiming(cds, GTOstim, 0);
[sequenceTimes, eventTimes] = getSequenceTimes(cds, stimState,GTOstim);
% use experiment class to generate PSTH
preTime = eventTimes(1) - sequenceTimes(1,1);
postTime = sequenceTimes(1,2)-eventTimes(1);

% stack events and spikes
spikes = [];
for i = 1:length(eventTimes)
    spikeMask = (cds.units(neuronNumber).spikes.ts  > eventTimes(i) - preTime & ...
        cds.units(neuronNumber).spikes.ts < eventTimes(i) + postTime);
    spikes = [spikes; cds.units(neuronNumber).spikes.ts(spikeMask) - eventTimes(i)];
end

% initialize bins
numBinsPre = ceil(preTime/userBinSize);
numBinsPost = ceil(postTime/userBinSize);
% set bin edges
binEdges = [0];
for i = 1:numBinsPre
    binEdges = [-i*userBinSize binEdges];
end
for i = 1:numBinsPost
    binEdges = [binEdges i*userBinSize];
end
[binCounts, binEdges] = histcounts(spikes,binEdges);
if(useRate)
    binCounts = binCounts/numel(eventTimes)/(binEdges(2)-binEdges(1));
    yLabelStr = 'Average Firing Rate (spikes/s)';
else
    binCounts = binCounts/numel(eventTimes);
    yLabelStr = 'Average Spike Counts (spikes)';
end
figure();
bar(binEdges(1:end-1) + mode(diff(binEdges))/2, binCounts);
% clean up graph
xlabel('Time After Stim Start (s)');
ylabel(yLabelStr);

% find optimal bin width if optimalBinSize is 1
if(optimalBinSize)
    figure();
    [optBinCounts, optBinEdges] = computeOptimalBinSize(spikes, userBinSize, preTime, postTime);
    if(useRate)
        optBinCounts = optBinCounts/numel(eventTimes)/(optBinEdges(2)-optBinEdges(1));
        yLabelStr = 'Average Firing Rate (spikes/s)';
    else
        optBinCounts = optBinCounts/numel(eventTimes);
        yLabelStr = 'Average Spike Counts (spikes)';
    end
    bar(optBinEdges(1:end-1) + + mode(diff(optBinEdges))/2, optBinCounts);
    xlabel('Time After Stim Start (s)');
    ylabel(yLabelStr);
end



% figure();
% [histogramOut] = PESTH(exp.units,eventTimes,preTime,postTime,neuronNumber);
% cd(funcFolder);
% warning('on')
% 
% plot([0.7 0.7],[0,max(get(gca,'YTick'))],'r')
end

% 
% exp = experiment();
% exp.meta.hasLfp=cds.meta.hasLfp;
% exp.meta.hasKinematics=cds.meta.hasKinematics;
% exp.meta.hasForce=cds.meta.hasForce;
% exp.meta.hasUnits=cds.meta.hasUnits;
% exp.meta.hasTrials=false;
% exp.meta.hasAnalog=cds.meta.hasAnalog;
% exp.addSession(cds);

