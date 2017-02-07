function [ barFig ] = plotPSTH(cds, neuronNumber, eventTimes, preTime, postTime, binSize, varargin)
% plots a PSTH given the data and relevant times/binSize. This should be
% usuable with any function that aims to plot a PSTH
useRate = 0;
gaussianSmooth = 0;
gaussianSmooth_std = binSize;
plotGaussianSmooth = 0;
eventOccurs = 0;
optimalBinSize = 0;
useEndAsZero = 0;
zeroCenter = 0;
averageSpikeWaveform = 0;

for i = 1:2:length(varargin)
    switch varargin{i}
        case 'optimalBinSize'
            optimalBinSize = varargin{i+1};
        case 'useRate'
            useRate = varargin{i+1};
        case 'eventOccurs'
            eventOccurs = varargin{i+1};
        case 'gaussianSmooth'
            gaussianSmooth = varargin{i+1};
        case 'gaussianStd'
            gaussianSmooth_std = varargin{i+1};
        case 'plotGaussian'
            plotGaussianSmooth = varargin{i+1};
        case 'useEndAsZero'
            useEndAsZero = varargin{i+1};
        case 'zeroCenter'
            zeroCenter = varargin{i+1};
        case 'averageSpikeWaveform'
            averageSpikeWaveform = varargin{i+1};
    end
    
end

% stack events and spikes
spikes = [];
for i = 1:length(eventTimes)
    spikeMask = (cds.units(neuronNumber).spikes.ts  > eventTimes(i) - preTime & ...
        cds.units(neuronNumber).spikes.ts < eventTimes(i) + postTime);
    spikes = [spikes; cds.units(neuronNumber).spikes.ts(spikeMask) - eventTimes(i)];
end

% find optimal bin width if optimalBinSize is 1
if(optimalBinSize)
    [optBinCounts, optBinEdges] = computeOptimalBinSize(spikes, binSize, preTime, postTime);
    binSize = optBinEdges(2) - optBinEdges(1);
end

% initialize bins and set edges
numBinsPre = ceil(preTime/binSize);
numBinsPost = ceil(postTime/binSize);
if(~zeroCenter)
    binEdges = [0];
    for i = 1:numBinsPre
        binEdges = [-i*binSize binEdges];
    end
    for i = 1:numBinsPost
        binEdges = [binEdges i*binSize];
    end
else
    binEdges = [-binSize/2 binSize/2];
    for i = 1:numBinsPre-1
        binEdges = [-i*binSize-binSize/2 binEdges];
    end
    for i = 1:numBinsPost-1
        binEdges = [binEdges i*binSize+binSize/2];
    end
end

% get bin counts and yLabelStr
if(eventOccurs)
    for i = 1:length(eventTimes)
        [eventInBin, binEdges] = histcounts(cds.units(neuronNumber).spikes.ts(spikeMask) - eventTimes(i),binEdges);
        if(i==1)
            binCounts = (eventInBin>0);
        else
            binCounts = binCounts+(eventInBin>0);
        end
    end
    binCounts = binCounts/length(eventTimes(i));
    yLabelStr = 'Rate of Spike Occuring';
elseif(useRate)
    [binCounts, binEdges] = histcounts(spikes,binEdges);
    binCounts = binCounts/numel(eventTimes)/(binEdges(2)-binEdges(1));
    yLabelStr = 'Average Firing Rate (spikes/s)';
else
    [binCounts, binEdges] = histcounts(spikes,binEdges);
    binCounts = binCounts/numel(eventTimes);
    yLabelStr = 'Average Spike Counts (spikes)';
end


% perform gaussian smoothing if required
if(gaussianSmooth)
    kernel_width = ceil(3*gaussianSmooth_std/binSize);
    kernel = normpdf(-kernel_width*binSize: ...
        binSize: ...
        kernel_width*binSize,...
        0, gaussianSmooth_std); 
    normalizer = conv(kernel,ones(1,length(binEdges)-1));
    smoothed_fr_inter = conv(kernel,binCounts)./normalizer;
    smoothed_fr = smoothed_fr_inter(kernel_width+1:end-kernel_width);
end


barFig = figure();
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

if(averageSpikeWaveform)
    % want 50-100 ms bin, this needs to be changed in the code
    binIdx = numBinsPre + 2;
    for i = 1:length(eventTimes)
        spikemask = (cds.units(neuronNumber).spikes.ts  > eventTimes(i) + binEdges(binIdx) & ...
            cds.units(neuronNumber).spikes.ts < eventTimes(i) + binEdges(binIdx+1));
        spikeIdx = find(cds.units(neuronNumber).spikes.ts(spikeMask));
    end
    spikesTable = cds.units(neuronNumber).spikes;
    waves = spikesTable{spikeIdx,:};
    aveWaves = mean(waves);
    figure();
    plot(aveWaves,'r','linewidth',2);
    hold on
    plot(aveWaves + std(waves),'--r','linewidth',2);
    plot(aveWaves - std(waves),'--r','linewidth',2);
    ylabel('Voltage (\muV)');
    xlabel('Wave Point');
end

end