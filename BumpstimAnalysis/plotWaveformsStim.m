function [  ] = plotWaveformsStim( cds, neuronNumber,chanNum,figNum,varargin )
maxWavesPlot = 100;

timeAfterStimRawNoStim = 20/1000;
timeAfterStimRawArtifact = 5/1000;
timeBeforeStimRawArtifact = 0/1000;
makeFigure = 1;
plotTitle = 0;
alignWaves = 1;
titleToPlot = num2str(neuronNumber);
waveformsSentExist = any(isfield(cds,'waveforms'));
stimElectrode = -1;
numChans = 1;
adjustForResolution = 1;

% save stuff
saveFigures = 0;
figDir = '';
figPrefix = '';

for i = 1:2:size(varargin,2)
    switch varargin{i}
        case 'timeAfterStimRawNoStim'
            timeAfterStimRawNoStim = varargin{i+1};
        case 'timeAfterStimRawArtifact'
            timeAfterStimRawArtifact = varargin{i+1};
        case 'timeBeforeStimRawArtifact'
            timeBeforeStimRawArtifact = varargin{i+1};
        case 'makeFigure'
            makeFigure = varargin{i+1};
        case 'plotTitle'
            plotTitle = varargin{i+1};
        case 'title'
            titleToPlot = varargin{i+1};
        case 'alignWaves'
            alignWaves = varargin{i+1};
        case 'stimElectrode'
            stimElectrode = varargin{i+1};
        case 'adjustForResolution'
            adjustForResolution = varargin{i+1};
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

% deals with making figure -- no fancy subplot stuff here
if(makeFigure)
    figure();
end

% chan list if necessary
if(any(isfield(cds.waveforms,'chanSent')))
    chanList = unique(cds.waveforms.chanSent);
    numChans = numel(chanList);
else
    chanList = stimElectrode;
    numChans = 1;
end

% get and plot actual data now
subplot(2,1,1) % raw waves not near artifact

wavesPlot = [];
spikeMask = zeros(numel(cds.units(neuronNumber).spikes.ts,1));
for st = 1:numel(cds.stimOn)
    if((~waveformsSentExist || cds.waveforms.waveSent(st)==figNum) && (~any(isfield(cds.waveforms,'chanSent')) || cds.waveforms.chanSent(st)==chanList(chanNum)))
        if(st==numel(cds.stimOn))
            spikeMask = spikeMask | (cds.units(neuronNumber).spikes.ts > cds.stimOn(st) + timeAfterStimRawNoStim);
        else
            spikeMask = spikeMask | (cds.units(neuronNumber).spikes.ts > cds.stimOn(st) + timeAfterStimRawNoStim & cds.units(neuronNumber).spikes.ts < cds.stimOn(st+1));
        end
    end
end

waveIdx = find(spikeMask);
waveIdx = datasample(waveIdx,min(maxWavesPlot,numel(waveIdx)),'Replace',false);
if(numel(waveIdx) > 0)
    for wave = 1:numel(waveIdx)
        rawIdx = getRawDataIdx(cds.units(neuronNumber).spikes.ts(waveIdx(wave)),cds.units(neuronNumber).chan,cds.rawData.ts,cds.rawData.elec);
        if(rawIdx ~= -1)
            wavesPlot = [wavesPlot;cds.rawData.waveforms(rawIdx,:)];
        end
    end
    xDataWaves = ((1:size(wavesPlot,2))-1)/30; % in ms
    if(adjustForResolution)
        wavesPlot = wavesPlot/0.254;
    end
    plot(xDataWaves,(wavesPlot-mean(wavesPlot(:,40:end),2)))
end

ylim([-800 800])

% deals with title requests
% deals with title requests
if(plotTitle)
    if(strcmp(titleToPlot,'') == 0)
        title(titleToPlot);
    elseif(numChans > 1 && waveformsSentExist)
        title(strcat('Stim Chan: ',num2str(chanList(chanNum)),' Wave: ',num2str(figNum)));
    elseif(numChans == 1 && waveformsSentExist)
        title(strcat('Stim Chan: ',num2str(stimElectrode),' Wave: ',num2str(figNum)));
    else

    end
end

ylabel('Voltage (\muV)')
formatForLee(gcf);

subplot(2,1,2) % waveforms near stimulation artifact

wavesPlot = [];
spikeMask = zeros(numel(cds.units(neuronNumber).spikes.ts,1));
for st = 1:numel(cds.stimOn)
    if((~waveformsSentExist || cds.waveforms.waveSent(st)==figNum) && (~any(isfield(cds.waveforms,'chanSent')) || cds.waveforms.chanSent(st)==chanList(chanNum)))
        spikeMask = spikeMask | (cds.units(neuronNumber).spikes.ts < cds.stimOn(st) + timeAfterStimRawArtifact & cds.units(neuronNumber).spikes.ts > cds.stimOn(st) + timeBeforeStimRawArtifact);
    end
end

waveIdx = find(spikeMask);
waveIdx = datasample(waveIdx,min(maxWavesPlot,numel(waveIdx)),'Replace',false);
if(numel(waveIdx) > 0)
    for wave = 1:min(maxWavesPlot,numel(waveIdx))
        rawIdx = getRawDataIdx(cds.units(neuronNumber).spikes.ts(waveIdx(wave)),cds.units(neuronNumber).chan,cds.rawData.ts,cds.rawData.elec);
        if(rawIdx ~= -1)
            wavesPlot = [wavesPlot;cds.rawData.waveforms(rawIdx,:)];
        end
    end
    xDataWaves = ((1:size(wavesPlot,2))-1)/30; % in ms
    if(adjustForResolution)
        wavesPlot = wavesPlot/0.254;
    end
    if(size(wavesPlot,1) == size(wavesPlot,2))
        plot(xDataWaves,(wavesPlot-mean(wavesPlot(:,40:end),2))')
    else
        plot(xDataWaves,(wavesPlot-mean(wavesPlot(:,40:end),2)))
    end
end

ylim([-800 800])
ylabel('Voltage (\muV)')
xlabel('Time (ms)')
formatForLee(gcf);

if(saveFigures)
    fname = strcat(figPrefix,'nn',num2str(neuronNumber),'_chan',num2str(cds.units(neuronNumber).chan),'_stimChan',num2str(chanList(chanNum)),'_waveNum',num2str(figNum),'_rawWaveforms');
    saveFiguresLab(gcf,figDir,fname);
end

end

