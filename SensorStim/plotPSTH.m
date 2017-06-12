function [ binEdges, binCounts, stdDiffCounts ] = plotPSTH(cds, neuronNumber, sequenceTimes, eventTimes, preTime,postTime,stimState, binSize, varargin)
% plots a PSTH given the data and relevant times/binSize. This should be
% usuable with any function that aims to plot a PSTH
spindleStim = 0;
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
plotEndStimulation = 0;
binsAbove = 0;
highlightBin = 1;
figMake = 1;
legendMake = 1;
legendFontSize = 12;
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'spindleStim'
            spindleStim = varargin{i+1};
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
        case 'plotEndStimulation'
            plotEndStimulation = varargin{i+1};
        case 'binsAbove'
            binsAbove = varargin{i+1};
        case 'figMake'
            figMake = varargin{i+1};
        case 'legendMake'
            legendMake = varargin{i+1};
        case 'highlightBin'
            highlightBin = varargin{i+1};
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
binEdgesAbove = [];
binSizeAbove = 0.005; % 10 ms bin size always
binCare = numBinsPre + 2;
if(~zeroCenter)
    binEdges = [0];
    for i = 1:numBinsPre
        binEdges = [-i*binSize binEdges];
    end
    for i = 1:numBinsPost
        binEdges = [binEdges i*binSize];
    end
    if(binsAbove == 1)
        binEdgesAbove = binSize:binSizeAbove:(binSize*2);
    end
else
    binEdges = [-binSize/2 binSize/2];
    for i = 1:numBinsPre-1
        binEdges = [-i*binSize-binSize/2 binEdges];
    end
    for i = 1:numBinsPost-1
        binEdges = [binEdges i*binSize+binSize/2];
    end
    if(binsAbove == 1)
        binEdgesAbove = binSize/2:binSizeAbove:binSize*3/2;
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
if(binsAbove)
    [binCountsAbove,] = histcounts(spikes,binEdgesAbove); % bin spikes
    binCountsAbove = binCountsAbove/sum(binCountsAbove); % make into a percent of spikes
end
% % get stdDiffCounts
stdDiffCounts = 0;
diffCounts = zeros(numel(eventTimes),2);
for i = 1:length(eventTimes)
        [eventInBin, binEdges] = histcounts(cds.units(neuronNumber).spikes.ts(spikeMask) - eventTimes(i),binEdges);
        diffCounts(i,1) =  eventInBin(binCare) - mean(eventInBin(1:binCare-1));
        diffCounts(i,2) = eventInBin(binCare) - mean(eventInBin(binCare+1:end));
end
stdDiffCounts = std(diffCounts./(binSize));
for i = 1:2
        if(stdDiffCounts(1,i) == 0) % either good or bad
            if(sum(diffCounts(:,i) <= 0))
                stdDiffCounts(1,i) = 10000;
            end
        end
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
    if(figMake)
        barFig = figure();
    else
        barFig = gcf;
    end
    barWidth = 0.95;
    if(eventOccurs)
        bar(1000*(binEdges(1:end-1) + mode(diff(binEdges))/2),binEventRate,...
            'barWidth',barWidth,'facecolor','k','edgecolor','none');
        ylim([0,1]);
    else
        bar(1000*(binEdges(1:end-1) + mode(diff(binEdges))/2), binCounts,...
            'barWidth',barWidth,'facecolor','k','edgecolor','none');
    end
    xlim([binEdges(1)*1000, binEdges(end)*1000]);
    greyIdx = 1;
    if(highlightBin)
        % color bar we care about a different color
        aHand = barFig.CurrentAxes;
        zeroIdx = find(binEdges==0);
        xData = aHand.Children.XData;
        yData = aHand.Children.YData*0;
        yData(zeroIdx+1) = aHand.Children.YData(zeroIdx+1);
        hold on
        bf = bar(xData,yData,...
            'parent',aHand,'faceColor','r','EdgeColor','none','barWidth',barWidth);
        set(gca,'XLim',aHand.XLim,'YLim',aHand.YLim)
        greyIdx = 2;
    end
    
    % color bar right after zero idx a different color
    aHand = barFig.CurrentAxes;
    zeroIdx = find(binEdges==0);
    xData = aHand.Children(greyIdx).XData;
    yData = aHand.Children(greyIdx).YData*0;
    yData(zeroIdx) = aHand.Children(greyIdx).YData(zeroIdx);
    hold on
    bf = bar(xData,yData,...
        'parent',aHand,'faceColor',[0.25 0.25 0.25],'EdgeColor','none','barWidth',barWidth);
    set(gca,'XLim',aHand.XLim,'YLim',aHand.YLim)
    
    % clean up graph
    if(useEndAsZero)
        xlabel('Time After Stimulation Ends (ms)')
    else
        xlabel('Time After Stimulation Onset (ms)');
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
        
        if(spindleStim)
            [meanConf, plusMinus] = bootstrapConfidenceIntervalSpindle(cds, neuronNumber, sequenceTimes, eventTimes, binSize,stimState);
        else
            [meanConf, plusMinus] = bootstrapConfidenceInterval(cds, neuronNumber, sequenceTimes, eventTimes, binSize);
        end
        if(useRate)
            meanConf = meanConf/binSize;
            plusMinus = plusMinus/binSize;
        end
        h1=plot([-10000,10000],[plusMinus(2),plusMinus(2)],'--','color','g','linewidth',3)
        
%         if(useRate)
%             poissVal = poissinv(0.975,meanCounts*(binSize)*numel(eventTimes))/numel(eventTimes)/binSize;
%         else
%             poissVal = poissinv(0.975,meanCounts*numel(eventTimes))/numel(eventTimes);
%         end
%         h3=plot([-100,100],[poissVal,poissVal],'m','linewidth',2);
        
        a.XLim = xLimits;
        a.YLim(2) = ceil(max(a.YLim(2)*1.3,plusMinus(2)*1.3));
        if(legendMake && ~plotEndStimulation)
            l=legend([h1],'99% confidence interval');
            set(l,'box','off','FontSize',legendFontSize,'location','northwest');
        end
    end
    formatForLee(gcf);
    if(plotEndStimulation == 1) % plot end of stimulation marker (average)
        stimSampRate = 1;
        if(~isempty(cds.analog))
            stimSampRate = 1/(cds.analog{1,1}.t(2) - cds.analog{1,1}.t(1));
        else
            stimSampRate = 1/(cds.lfp.t(2) - cds.lfp.t(1));
        end
        stimStateTimeAverage = (sum(stimState == 1))/numel(eventTimes)/stimSampRate;
        
        % plot vertical line at this time
        hold on
        figHandle = gcf;
        hr = plot([stimStateTimeAverage*1000, stimStateTimeAverage*1000],figHandle.CurrentAxes.YLim,'--r','linewidth',3)
        if(legendMake && ~confInter)
            l=legend([hr],'Average Stimulation End');
            set(l,'box','off','FontSize',legendFontSize,'location','northwest');
        else
            [l,objh]=legend([h1,hr],'99% confidence interval','Average Stimulation End');
            objh(3).XData = [-0.1 0.1830];
            objh(5).XData = [-0.1 0.1830];
            set(l,'box','off','FontSize',legendFontSize,'location','northwest');
            p = l.Position;
            l.Position = p + [0.1 0 0 0];
        end
    end
    
    if(binsAbove)
        % using binCountsAbove and binEdgesAbove, plot a bar graph above
        % the bin we care about....this might suck to write
        % its probably easiest to just plot lines?
        
        % figure out line height based on max bin height
        a=gca;
        yl = a.YLim;
        maxLineHeight = 1/3*max(binCounts);
        % set yLim to compensate
        if(a.YLim(2) < (1/3+0.02)*max(binCounts) + max(binCounts))
            a.YLim(2) = (1/3+0.02)*max(binCounts) + 1.02*max(binCounts);
        end
        hold on
        
        % plot lines lol
        binEdgesAbove = binEdgesAbove(1:end-1) + (binEdgesAbove(2)-binEdgesAbove(1))/2;
        for i = 1:numel(binCountsAbove)
            lineHeight = maxLineHeight*binCountsAbove(i);
            xPos = [binEdgesAbove(i)*1000,1000*binEdgesAbove(i)];
            yPos = [max(binCounts)*1.02,max(binCounts)*1.02+lineHeight];
            plot(xPos,yPos,'r','lineWidth',1.5);         
        end
        
        a.YLim = yl;
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
        if(size(wavesWave,1)>1)
            aveWavesWave = mean(wavesWave);
        else
            aveWavesWave = wavesWave;
        end
        if(size(wavesAllElse > 1))
            aveWavesAllElse = mean(wavesAllElse);
        else
            aveWavesAllElse = wavesAllElse;
        end
        
        if(figMake)
            figure();
        end
        xWaves = (0:1:numel(aveWavesWave)-1)/30000*1000;
        plot(xWaves,aveWavesWave,'r','linewidth',2);
        hold on
        plot(xWaves,aveWavesAllElse,'k','linewidth',2);
        plot(xWaves,aveWavesWave + std(wavesWave),'--r','linewidth',2);
        plot(xWaves,aveWavesWave - std(wavesWave),'--r','linewidth',2);
        plot(xWaves,aveWavesAllElse + std(wavesAllElse),'--k','linewidth',2);
        plot(xWaves,aveWavesAllElse - std(wavesAllElse),'--k','linewidth',2);
        ylabel('Voltage (\muV)');
        xlabel('Time (ms)');
        l=legend('40-80ms waveforms','all other waveforms');
        set(l,'box','off','fontsize',legendFontSize,'location','best');
        formatForLee(gcf);
    end
end

end