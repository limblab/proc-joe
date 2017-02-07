function [ output_args ] = generatePSTHPeaks( cds, varargin )
% bins around the peaks of the spindle stimulation, does a gaussian
% smoothing if requested
%           'gaussianSmooth', # <- 0,1 

gaussianSmooth = 0;
GTOstim = 0;
useEndAsZero = 0;

for i = 1:2:length(varargin)
    switch varargin
        case 'gaussianSmooth'
            gaussianSmooth = varargin{i+1};
    end
end


% get stim timing
[stimState,] = determineStimTiming(cds, GTOstim, 0);
[spindlePeakTimes] = getSpindleStimPeaks(cds,stimState);
% [sequenceTimes, eventTimes] = getSequenceTimes(cds, stimState,GTOstim,useEndAsZero);


end

