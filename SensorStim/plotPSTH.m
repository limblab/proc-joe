function [ binEdges, binCounts ] = plotPSTH(cds, neuronNumber, eventTimes, preTime, postTime, binSize, varargin)
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
noPlots = 0;
confInter = 1;
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
        case 'noPlots'
            noPlots = varargin{i+1};
        case 'confidenceInterval'
            confInter = varargin{i+1};
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
    binCounts = binCounts/numel(eventTimes);
    yLabelStr = 'Rate of Spike Occuring';
elseif(useRate)
    [binCounts, binEdges] = histcounts(spikes,binEdges);
    binCounts = binCounts/numel(eventTimes)/(binSize);
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

if(~noPlots)
    barFig = figure();
    if(eventOccurs)
        bar(binEdges(1:end-1) + mode(diff(binEdges))/2,binEventRate);
        ylim([0,1]);
    elseif(useRate)
        bar(binEdges(1:end-1) + mode(diff(binEdges))/2, binCounts);
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
    
    %% conf interval if requested
    if(confInter) % so definitely not a normal distribution -- bootstrapping
        a = gca;
        xLimits = a.XLim;
        hold on
        [meanConf, plusMinus] = bootstrapConfidenceInterval(cds, neuronNumber, eventTimes(1), binSize);
        if(useRate)
            meanConf = meanConf/binSize;
            plusMinus = plusMinus/binSize;
        end
        plot([-100,100],[meanConf,meanConf],'k','linewidth',2)
        plot([-100,100],[plusMinus(1),plusMinus(1)],'k','linewidth',2)
        plot([-100,100],[plusMinus(2),plusMinus(2)],'k','linewidth',2)
        a.XLim = xLimits;
    end
    
    
    %% Average Spike Waveform
    spikeIdxWave = [];
    spikeIdxAllElse = [];
    
    temp = num2str(averageSpikeWaveform);
    
    if(str2num(temp(1))==1)
        % 
        if(length(temp)>1)
            binIdx = numBinsPre + str2num(temp(2:end));
        else
            binIdx = numBinsPre + 2;
        end
        for i = 1:length(eventTimes)
            spikeMask = (cds.units(neuronNumber).spikes.ts  > eventTimes(i) + binEdges(binIdx) & ...
                cds.units(neuronNumber).spikes.ts < eventTimes(i) + binEdges(binIdx+1));
            spikeIdxWave = [spikeIdxWave; find(spikeMask)];
            spikeIdxAllElse = [spikeIdxAllElse; find(~spikeMask)];
        end
        spikesTable = cds.units(neuronNumber).spikes;
        wavesWave = spikesTable{spikeIdxWave,:};
        wavesAllElse = spikesTable{spikeIdxAllElse,:};
        aveWavesWave = mean(wavesWave);
        aveWavesAllElse = mean(wavesAllElse);
        
        figure();
        plot(aveWavesWave,'r','linewidth',2);
        hold on
        plot(aveWavesAllElse,'b','linewidth',2);
        plot(aveWavesWave + std(wavesWave),'--r','linewidth',2);
        plot(aveWavesWave - std(wavesWave),'--r','linewidth',2);
        plot(aveWavesAllElse + std(wavesAllElse),'--b','linewidth',2);
        plot(aveWavesAllElse - std(wavesAllElse),'--b','linewidth',2);
        ylabel('Voltage (\muV)');
        xlabel('Wave Point');
        legend('binned waves','all other waves')
    end
end

end