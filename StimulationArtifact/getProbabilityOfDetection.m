function [probOut] = getProbabilityOfDetection(cds,neuronNumber,varargin)

lowTime = 0/1000;
highTime = 5/1000;
preTime = 20/1000;
postTime = 60/1000;

waveformsSentExist = any(isfield(cds,'waveforms'));
waveformTypes = 1;
if(waveformsSentExist)
    waveformTypes = 1:1:numel(unique(cds.waveforms.waveSent));
end
chans = 1;
stimElectrode = -1;
bumpTask = 0;
alignWaves = 1;
for i = 1:2:numel(varargin)
    switch varargin{i}
        case 'peakPeriod'
            temp = varargin{i+1};
            lowTime = temp(1);
            highTime = temp(2);
        case 'waveformTypes'
            waveformTypes = varargin{i+1};
        case 'chans'
            chans = varargin{i+1};
        case 'stimElectrode'
            stimElectrode = varargin{i+1};
        case 'preTime'
            preTime = varargin{i+1};
        case 'postTime'
            postTime = varargin{i+1};
    end
end

%% extract number of waveform types if applicable
if(waveformsSentExist)
    numWaveformTypes = numel(unique(cds.waveforms.waveSent));
end

if(any(isfield(cds.waveforms,'chanSent')))
    numChans = numel(unique(cds.waveforms.chanSent));
    chanList = unique(cds.waveforms.chanSent);
else
    chanList = stimElectrode;
    numChans = 1;
end

%% extract data
for c = 1:numChans
    for i = 1:numWaveformTypes
        spikeTimeData{c,i} = [];
    end
end
stimNum = zeros(numChans,numWaveformTypes);
if(bumpTask) % plot things aligned to bump times and whatnot
    % write this later :D
else % get data after stimulations
    for st = 1:numel(cds.stimOn)
        spikeMask = cds.units(neuronNumber).spikes.ts > cds.stimOn(st)-preTime & cds.units(neuronNumber).spikes.ts < cds.stimOn(st)+postTime;
        spikesPlot = (cds.units(neuronNumber).spikes.ts(spikeMask) - cds.stimOn(st));
        if(alignWaves)
            spikesPlot = spikesPlot + 4/30000;
        end
        numWaves = sum(spikeMask==1);
        if(waveformsSentExist)
            if(any(isfield(cds.waveforms,'chanSent')))
                chanNumber = find(unique(cds.waveforms.chanSent)==cds.waveforms.chanSent(st));
                stimNum(chanNumber,cds.waveforms.waveSent(st)) = stimNum(chanNumber,cds.waveforms.waveSent(st)) + 1;
                if(~isempty(spikesPlot))
                    spikeTimeData{chanNumber,cds.waveforms.waveSent(st)}(end+1:end+numWaves) = spikesPlot';
                end
            else
                stimNum(1,cds.waveforms.waveSent(st)) = stimNum(1,cds.waveforms.waveSent(st)) + 1;
                if(~isempty(spikesPlot))
                    spikeTimeData{1,cds.waveforms.waveSent(st)}(end+1:end+numWaves) = spikesPlot';
                end
            end
        else
            if(any(isfield(cds.waveforms,'chanSent')))
                chanNumber = find(unique(cds.waveforms.chanSent)==cds.waveforms.chanSent(st));
                stimNum(chanNumber,1) = stimNum(chanNumber,1) + 1;
                spikeTimeData{chanNumber,1}(end+1:end+numWaves) = spikesPlot';
            else
                stimNum(1,1) = stimNum(1,1) + 1;
                spikeTimeData{1,1}(end+1:end+numWaves) = spikesPlot';
            end
        end
    end
end

%% calculate probability of detection
prob = zeros(numChans,numWaveformTypes);
for c = 1:numChans
    for w = 1:numWaveformTypes
        baseline = sum(spikeTimeData{c,w} < min(0,lowTime))/(min(0,lowTime)+preTime); % spikes/ms
        ROI = sum(spikeTimeData{c,w} > lowTime & spikeTimeData{c,w} < highTime)/(highTime-lowTime); % spikes/ms
        prob(c,w) = (ROI - baseline)/numel(cds.stimOn);
    end
end

probOut = prob(chans,waveformTypes);

end