function [ histogramOut ] = generatePESTH(cds,neuronNumber, varargin)
% possible varargin: 'useRate', 1 <- plots rate 
%                    'optimalBinSize', 1 <-calculates and plots PSTH with optimal bin size as well, 
%                    'binSize', # <- uses this binSize instead of 0.05.
%                    'zeroEvent' 'start','end' <- zero point of histogram
%                    'gaussianSmooth', <- 1 if yes, 0 if no
%                    'gaussianStd', # <- # is the std of the kernel
%                    'plotGaussianSmooth', # <- 1 or 0 if want to plot
% initialize variables
userBinSize = 0.05;
optimalBinSize = 0;
useRate = 0;
useEndAsZero = 0;
gaussianSmooth = 0;
gaussianSmooth_std = userBinSize;
plotGaussianSmooth = 0;
eventOccurs = 0;
% deal with varargin
for i = 1:2:size(varargin,2)
    switch varargin{i}
        case 'useRate'
            useRate = varargin{i+1};
        case 'optimalBinSize'
            optimalBinSize = varargin{i+1};
        case 'binSize'
            userBinSize = varargin{i+1};
        case 'zeroEvent'
            if(strcmp(lower(varargin{i+1}), 'end'))
                useEndAsZero = 1;
            end
        case 'gaussianSmooth'
            gaussianSmooth = varargin{i+1};
        case 'gaussianStd'
            gaussianSmooth_std = varargin{i+1};
        case 'plotGaussian'
            plotGaussianSmooth = varargin{i+1};
        case 'eventOccurs'
            eventOccurs = varargin{i+1};
    end
end

% generates a post stimulus time histogram from the given unit in the cds
GTOstim = 0;
if(size(cds.analog,1) == 0)
    GTOstim = 1;
end

[stimState,] = determineStimTiming(cds, GTOstim, 0);
[sequenceTimes, eventTimes] = getSequenceTimes(cds, stimState,GTOstim,useEndAsZero);
% use experiment class to generate PSTH
preTime = eventTimes(1) - sequenceTimes(1,1);
postTime = sequenceTimes(1,2)-eventTimes(1);

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

% stack events and spikes
spikes = [];
binEventOccurences = [];
for i = 1:length(eventTimes)
    spikeMask = (cds.units(neuronNumber).spikes.ts  > eventTimes(i) - preTime & ...
        cds.units(neuronNumber).spikes.ts < eventTimes(i) + postTime);
    spikes = [spikes; cds.units(neuronNumber).spikes.ts(spikeMask) - eventTimes(i)];
    % count times at least one spike happens per bin
    [eventInBin, binEdges] = histcounts(cds.units(neuronNumber).spikes.ts(spikeMask) - eventTimes(i),binEdges);
    if(i==1)
        binEventOccurences = (eventInBin>0);
    else
        binEventOccurences = binEventOccurences+(eventInBin>0);
    end
end
binEventRate = binEventOccurences./length(eventTimes);


[binCounts, binEdges] = histcounts(spikes,binEdges);
if(eventOccurs)
    yLabelStr = 'Rate of Spike Occuring';
elseif(useRate)
    binCounts = binCounts/numel(eventTimes)/(binEdges(2)-binEdges(1));
    yLabelStr = 'Average Firing Rate (spikes/s)';
else
    binCounts = binCounts/numel(eventTimes);
    yLabelStr = 'Average Spike Counts (spikes)';
end

% perform gaussian smoothing if required
if(gaussianSmooth)
    kernel_width = ceil(3*gaussianSmooth_std/userBinSize);
    kernel = normpdf(-kernel_width*userBinSize: ...
        userBinSize: ...
        kernel_width*userBinSize,...
        0, gaussianSmooth_std); 
    normalizer = conv(kernel,ones(1,length(binEdges)-1));
    smoothed_fr_inter = conv(kernel,binCounts)./normalizer;
    smoothed_fr = smoothed_fr_inter(kernel_width+1:end-kernel_width);
end


figure();
if(eventOccurs)
    bar(binEdges(1:end-1) + mode(diff(binEdges))/2,binEventRate);
    ylim([0,1]);
elseif(useRate)
    bar(binEdges(1:end-1) + mode(diff(binEdges))/2, binCounts);
    ylim([0,1]);
else
    bar(binEdges(1:end-1) + mode(diff(binEdges))/2, binCounts);
end
% clean up graph
if(useEndAsZero)
    xlabel('Time Relative to Stim End (s)')
else
    xlabel('Time Relative to Stim Start (s)');
end
ylabel(yLabelStr);
if(plotGaussianSmooth && gaussianSmooth)
    hold on
    plot((binEdges(1:end-1)+binEdges(2:end))/2,smoothed_fr,'r','linewidth',2);
end
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

end

