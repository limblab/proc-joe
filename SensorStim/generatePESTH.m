function [binEdges,binCounts] = generatePESTH(cds,neuronNumber, varargin)
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
        case 'confidenceInterval'
            confidenceInterval = varargin{i+1};
    end
end

% generates a post stimulus time histogram from the given unit in the cds
GTOstim = 0;
if(size(cds.analog,1) == 0)
    GTOstim = 1;
end

if(~exist('sequenceTimes') || ~exist('eventTimes'))
    [stimState,] = determineStimTiming(cds, GTOstim, 0);
    [sequenceTimes, eventTimes] = getSequenceTimes(cds, stimState,GTOstim,useEndAsZero);
end
preTime = eventTimes(1) - sequenceTimes(1,1);
postTime = sequenceTimes(1,2)-eventTimes(1);



[binEdges, binCounts] = plotPSTH(cds, neuronNumber, eventTimes, preTime, postTime, userBinSize,...
    'optimalBinSize', optimalBinSize, 'useRate',useRate,...
    'eventOccurs',eventOccurs,'gaussianSmooth',gaussianSmooth,...
    'gaussianStd',gaussianSmooth_std,'plotGaussian',plotGaussianSmooth,...
    'useEndAsZero',useEndAsZero,'averageSpikeWaveform',averageSpikeWaveform,...
    'plotWaveform',plotWaveform, 'noPlots',noPlots);
end

