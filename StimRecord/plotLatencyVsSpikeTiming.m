function [] = plotLatencyVsSpikeTiming(cds,neuronNumber,varargin)

timeAfterStim = 5/1000;

plotTitle = 1;
titleToPlot = '';
makeFigure = 1;
waveformsMakeSubplots = 0;

waveformsSentExist = any(isfield(cds,'waveforms'));
waveformsMakeSubplots = 0;
waveformTypesPlot = 1;
if(waveformsSentExist)
    waveformTypesPlot = 1:1:numel(unique(cds.waveforms.waveSent));
end

%% deal with varagin
for i = 1:2:size(varargin,2)
    switch varargin{i}
        case 'timeAfterStim'
            timeAfterStim = varargin{i+1};
        case 'plotTitle'
            plotTitle = varargin{i+1};
        case 'title'
            titleToPlot = varargin{i+1};
        case 'makeFigure'
            makeFigure = varargin{i+1};
        case 'makeSubplots'
            waveformsMakeSubplots = varargin{i+1};
        case 'waveformTypes'
            waveformTypesPlot = varargin{i+1};
    end
end

%% extract number of waveform types if applicable
if(waveformsSentExist)
    numWaveformTypes = numel(unique(cds.waveforms.waveSent));
end

%% grab and plot data

for fig = waveformTypesPlot
    if(makeFigure)
        if(waveformsSentExist && waveformsMakeSubplots && fig == 1)
            figure(); % make figure for the subplots
            subplot(numWaveformTypes,1,fig);
        elseif(waveformsSentExist && waveformsMakeSubplots)
            subplot(numWaveformTypes,1,fig); % subplots
        elseif(waveformsSentExist && ~waveformsMakeSubplots)
            figure(); % figure for each waveform type
        elseif(fig==1)
            figure(); % figure for a single waveform type
        end
    end
    
    latency = [];
    timeBetween = [];
    for st = 1:numel(cds.stimOn)
        if(~waveformsSentExist || (waveformsSentExist && cds.waveforms.waveSent(st) == fig))
            spikeMask = (cds.units(neuronNumber).spikes.ts < cds.stimOn(st) + timeAfterStim & cds.units(neuronNumber).spikes.ts > cds.stimOn(st));
            spikeIdx = find(spikeMask);
            if(numel(spikeIdx) > 0)
                spikeIdx = spikeIdx(1); % get first one
                if(spikeIdx ~= 1)
                    spikeTs = cds.units(neuronNumber).spikes.ts(spikeIdx) - cds.stimOn(st);
                    prevSpikeTs = cds.units(neuronNumber).spikes.ts(spikeIdx-1)-cds.stimOn(st);
                    latency(end+1:end+numel(spikeTs)) = spikeTs*1000;
                    timeBetween(end+1:end+numel(prevSpikeTs)) = (spikeTs - prevSpikeTs)*1000;
                end
            end
        end
    end


    plot(timeBetween,latency,'k.')
    ylabel('Time since last action potential (ms)')
    xlabel('Latency (ms)')
    % deals with title requests
    if(plotTitle)
        if(strcmp(titleToPlot,'') == 0)
            title(titleToPlot);
        elseif(waveformsSentExist)
            title(num2str(fig));
        else
            title(num2str(neuronNumber));
        end
    end
    formatForLee(gcf);
    
end

end