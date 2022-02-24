function [barFig] = generatePSTHPeaks( cds,neuronNumber, varargin )
% bins around the peaks of the spindle stimulation, does a gaussian
% smoothing if requested
%           'peak', + or - # for which peak to center around
%           'zeroCenter', 1 if bin centers around zero, 0 if zero is edge
gaussianSmooth = 0;
gaussianSmooth_std = 0.01;
GTOstim = 0;
useRate = 0;
plotGaussian = 0;
binSize = -1;
peak = 1;
zeroCenter = 0;
plotWaveform = 0;

for i = 1:2:length(varargin)
    switch varargin{i}
        case 'gaussianSmooth'
            gaussianSmooth = varargin{i+1};
        case 'gaussianStd'
            gaussianSmooth_std = varargin{i+1};
        case 'plotGaussian'
            plotGaussian = varargin{i+1};
        case 'binSize'
            binSize = varargin{i+1};
        case 'useRate'
            useRate = varargin{i+1};
        case 'peak'
            peak = sign(varargin{i+1});
        case 'zeroCenter'
            zeroCenter = varargin{i+1};
        case 'plotWaveform'
            plotWaveform = varargin{i+1};
    end
end


% get stim timing and peak times
[stimState,] = determineStimTiming(cds, GTOstim, 0);
[spindlePeakTimes,spindleStimTimes] = getSpindleStimPeaks(cds,stimState,peak);
frequencyStim = 1/(spindlePeakTimes(2)-spindlePeakTimes(1));
if(binSize == -1)
    binSize = (1/frequencyStim)/15;
end
% bin around peaks
preTime = (spindlePeakTimes(2)-spindlePeakTimes(1))/2;
postTime = preTime;

barFig = plotPSTH(cds, neuronNumber, spindlePeakTimes, preTime, postTime, binSize,...
    'useRate',useRate,'gaussianSmooth',gaussianSmooth,...
    'gaussianStd',gaussianSmooth_std,'plotGaussian',plotGaussian,...
    'zeroCenter',zeroCenter);
xlabel('Time (s)')
if(plotWaveform) % plots the stimulation waveform on top of bar graph
    % plot times = preTime to postTime, centered around a spindlePeakTime
    spindlePeakIndex = [find(cds.analog{1,1}.t==spindlePeakTimes(3)), find(cds.analog{1,1}.t==spindlePeakTimes(4))];
    plotIndex = [spindlePeakIndex(1)-ceil((spindlePeakIndex(2)-spindlePeakIndex(1))/2), ...
        ceil((spindlePeakIndex(1)+spindlePeakIndex(2))/2)];
    
    plotTimes = linspace(-preTime,postTime,plotIndex(2)-plotIndex(1)+1);
    data = cds.analog{1,1};
    plotVals = data{plotIndex(1):plotIndex(2),2};
    % plotVals normalized to 0 to magnitude/2 of bargraph
    plotVals = (plotVals-min(plotVals))./(max(plotVals)-min(plotVals))*max(barFig.Children.YLim);
    % shift plotVals above bargraph
%     plotVals = plotVals + max(barFig.Children.YLim)/2;
    hold on
    plot(plotTimes,plotVals,'-r','linewidth',3)
end

end

