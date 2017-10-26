function [  ] = plotInterspikeIntervalHistogram( cds, neuronNumber, varargin )
plotTitle = 0;
title = '';
makeFigure = 1;
binSize = 0.1; % in ms
xLimits = [0,10]; % in ms
displayText = 1;

%
stimElectrode = -1;

% save stuff
saveFigures = 0;
figDir = '';
figPrefix = '';

%% deal with varagin
for i = 1:2:size(varargin,2)
    switch varargin{i}
        case 'plotTitle'
            plotTitle = varargin{i+1};
        case 'title'
            titleToPlot = varargin{i+1};
        case 'makeFigure'
            makeFigure = varargin{i+1};
        case 'binSize'
            binSize = varargin{i+1};
        case 'xLim'
            xLimits = varargin{i+1};
        case 'displayText'
            displayText = varargin{i+1};
        case 'stimElectrode'
            stimElectrode = varargin{i+1};
        case 'saveFigures'
            saveFigures = varargin{i+1};
        case 'figDir'
            figDir = varargin{i+1};
        case 'figPrefix'
            figPrefix = varargin{i+1};
        
    end
end

if(saveFigures && strcmp(figDir,''))
    saveFigures = 0;
end

if(makeFigure)
    figure();
end

%% get data
interSpikeInterval = diff(cds.units(neuronNumber).spikes.ts);
interSpikeInterval = interSpikeInterval*1000; % seconds to ms

%% bin and plot data
binEdges = 0:binSize:max(xLimits);
[binCounts,~] = histcounts(interSpikeInterval,binEdges);
bar(binEdges(1:end-1)+(binEdges(2)-binEdges(1))/2,binCounts)
ylabel('Spike count')
xlabel('Time (ms)')
formatForLee(gcf)
%% display things (namely percentage below 1.7ms)
if(displayText)
    percentageBelow = sum(interSpikeInterval < 1.7)/numel(interSpikeInterval)*100;
    disp(percentageBelow)
end

if(saveFigures)
    fname = strcat(figPrefix,'nn',num2str(neuronNumber),'_chan',num2str(cds.units(neuronNumber).chan),'_ISI');
    saveFiguresLIB(gcf,figDir,fname);
end

end

