function [probActual,probExpected] = getDoublePeakProbability(cds,neuronNumber,varargin)

peakPeriod1 = [0,4];
peakPeriod2 = [5,9];
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

for i = 1:2:numel(varargin)
    switch varargin{i}
        case 'peakPeriod1'
            peakPeriod1 = varargin{i+1};
        case 'peakPeriod2'
            peakPeriod2 = varargin{i+1};
        case 'preTime'
            preTime = varargin{i+1};
        case 'postTime'
            postTime = varargin{i+1};
        case 'chans'
            chans = varargin{i+1};
        case 'waveformTypes'
            waveformTypes = varargin{i+1};
        case 'stimElectrode'
            stimElectrode = varargin{i+1};
        case 'bumpTask'
            bumpTask = varargin{i+1};
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
        spikePresentPeriod1{c,i} = [];
        spikePresentPeriod2{c,i} = [];
        spikePresentBoth{c,i} = [];
    end
end
stimNum = zeros(numChans,numWaveformTypes);
if(bumpTask) % plot things aligned to bump times and whatnot
    % write this later :D
else % get data after stimulations
    for st = 1:numel(cds.stimOn)
        % see if spike is present in first peak period
        spikeMask1 = cds.units(neuronNumber).spikes.ts > cds.stimOn(st)+peakPeriod1(1) & cds.units(neuronNumber).spikes.ts < cds.stimOn(st)+peakPeriod1(2);
        spikeMask2 = cds.units(neuronNumber).spikes.ts > cds.stimOn(st)+peakPeriod2(1) & cds.units(neuronNumber).spikes.ts < cds.stimOn(st)+peakPeriod2(2);

        if(waveformsSentExist)
            if(any(isfield(cds.waveforms,'chanSent')))
                chanNumber = find(unique(cds.waveforms.chanSent)==cds.waveforms.chanSent(st));
                stimNum(chanNumber,cds.waveforms.waveSent(st)) = stimNum(chanNumber,cds.waveforms.waveSent(st)) + 1;
                spikePresentPeriod1{chanNumber,cds.waveforms.waveSent(st)}(end+1,1) = sum(spikeMask1);
                spikePresentPeriod2{chanNumber,cds.waveforms.waveSent(st)}(end+1,1) = sum(spikeMask2);
                spikePresentBoth{chanNumber,cds.waveforms.waveSent(st)}(end+1,1) = sum(spikeMask1~=0) && sum(spikeMask2)~=0;
            else
                stimNum(1,cds.waveforms.waveSent(st)) = stimNum(1,cds.waveforms.waveSent(st)) + 1;
                spikePresentPeriod1{1,cds.waveforms.waveSent(st)}(end+1,1) = sum(spikeMask1)~=0;
                spikePresentPeriod2{1,cds.waveforms.waveSent(st)}(end+1,1) = sum(spikeMask2)~=0;
                spikePresentBoth{1,cds.waveforms.waveSent(st)}(end+1,1) = sum(spikeMask1~=0) && sum(spikeMask2)~=0;
            end
        else
            if(any(isfield(cds.waveforms,'chanSent')))
                chanNumber = find(unique(cds.waveforms.chanSent)==cds.waveforms.chanSent(st));
                stimNum(chanNumber,1) = stimNum(chanNumber,1) + 1;
                spikePresentPeriod1{chanNumber,1}(end+1,1) = sum(spikeMask1)~=0;
                spikePresentPeriod2{chanNumber,1}(end+1,1) = sum(spikeMask2)~=0;
                spikePresentBoth{chanNumber,1}(end+1,1) = sum(spikeMask1~=0) && sum(spikeMask2)~=0;
            else
                stimNum(1,1) = stimNum(1,1) + 1;
                spikePresentPeriod1{1,1}(end+1,1) = sum(spikeMask1)~=0;
                spikePresentPeriod2{1,1}(end+1,1) = sum(spikeMask2)~=0;
                spikePresentBoth{1,1}(end+1,1) = sum(spikeMask1~=0) && sum(spikeMask2)~=0;
            end
        end
    end
end


%% output probability of second spike given a first spike has happened
probActual = zeros(numChans,numWaveformTypes);
probExpected = zeros(numChans,numWaveformTypes);
for c = 1:numChans
    for w = 1:numWaveformTypes
        probActual(c,w) = sum(spikePresentBoth{c,w})/size(spikePresentBoth{c,w},1);
        probExpected(c,w) = sum(spikePresentPeriod1{c,w})/size(spikePresentPeriod1{c,w},1)*sum(spikePresentPeriod2{c,w})/size(spikePresentPeriod2{c,w},1);
    end
end
probActual = probActual(chans,waveformTypes);
probExpected = probExpected(chans,waveformTypes);

end