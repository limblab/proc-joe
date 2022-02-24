function [binEdges,binCounts,stdDiffCounts] = generatePESTH(cds,neuronNumber, varargin)
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
useRate = 1;
useEndAsZero = 0;
gaussianSmooth = 0;
gaussianSmooth_std = userBinSize;
plotGaussianSmooth = 0;
eventOccurs = 0;
averageSpikeWaveform = 0;
plotWaveform = 0;
noPlots=0;
confidenceInterval = 1;
plotEndStimulation = 0;
binsAbove = 0;
figMake = 1;
legendMake = 1;
highlightBin = 1;
spindleStim = 0;
% deal with varargin
for i = 1:2:size(varargin,2)
    switch varargin{i}
        case 'spindleStim'
            spindleStim = varargin{i+1};
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
            plotGaussianSmooth = gaussianSmooth;
        case 'gaussianStd'
            gaussianSmooth_std = varargin{i+1};
        case 'eventOccurs'
            eventOccurs = varargin{i+1};
        case 'averageSpikeWaveform'
            averageSpikeWaveform = varargin{i+1};
        case 'plotWaveform'
            plotWaveform = varargin{i+1};
        case 'noPlots'
            noPlots = varargin{i+1};
        case 'sequenceTimes'
            sequenceTimes = varargin{i+1};
        case 'eventTimes'
            eventTimes = varargin{i+1};
        case 'stimState'
            stimState = varargin{i+1};
        case 'confidenceInterval'
            confidenceInterval = varargin{i+1};
        case 'plotEndStimulation',
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

% generates a post stimulus time histogram from the given unit in the cds
GTOstim = 0;
if(size(cds.analog,1) == 0)
    GTOstim = 1;
end

if(~exist('sequenceTimes') || ~exist('eventTimes') || ~exist('stimState'))
    [stimState,] = determineStimTiming(cds, GTOstim, 0);
    [sequenceTimes, eventTimes] = getSequenceTimes(cds, stimState,GTOstim,useEndAsZero);
end
preTime = max(eventTimes(:) - sequenceTimes(:,1));
postTime = max(sequenceTimes(:,2)-eventTimes(:));



[binEdges, binCounts, stdDiffCounts] = plotPSTH(cds, neuronNumber, sequenceTimes, eventTimes, preTime, postTime, stimState, userBinSize,...
    'optimalBinSize', optimalBinSize, 'useRate',useRate,...
    'eventOccurs',eventOccurs,'gaussianSmooth',gaussianSmooth,...
    'gaussianStd',gaussianSmooth_std,'plotGaussian',plotGaussianSmooth,...
    'useEndAsZero',useEndAsZero,'averageSpikeWaveform',averageSpikeWaveform,...
    'plotWaveform',plotWaveform, 'noPlots',noPlots,'confidenceInterval',confidenceInterval,...
    'plotEndStimulation',plotEndStimulation,'binsAbove',binsAbove,'figMake',figMake,'legendMake',legendMake,...
    'highlightBin',highlightBin,'spindleStim',spindleStim);

end

