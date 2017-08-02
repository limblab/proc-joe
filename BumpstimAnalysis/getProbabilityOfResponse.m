function [probOut] = getProbabilityOfResponse(cds,neuronNumber,varargin)

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
autoPeakPeriod = 0;
for i = 1:2:numel(varargin)
    switch varargin{i}
        case 'peakPeriod'
            if(~strcmp('automatic',varargin{i+1}))
                temp = varargin{i+1};
                lowTime = temp(1);
                highTime = temp(2);
            else
                autoPeakPeriod = 1;
            end
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
        if(autoPeakPeriod) % get high time and low time automatically
            % bin data
            binSize = 0.0002;
            bE = -preTime:binSize:postTime;
            [bC,~] = histcounts(spikeTimeData{c,w},bE);
            % find peak in PSTH between 0 and 5ms
            zeroBin = floor(preTime/binSize);
            [~,peakIdx] = max(bC(zeroBin:zeroBin+ceil(0.005/binSize)));
            peakIdx = peakIdx + zeroBin - 1;
            % get baseline rate
            baselineRate = bC(1:zeroBin-floor(0.003/binSize)); % preTime to -3ms
            % set high and low time based on baseline firing rate
            lowTemp = find(bC(1:peakIdx) < baselineRate);
            lowTime = bE(lowTemp(end));
            highTemp = find(bC(peakIdx:end) < baselineRate);
            highTime = highTemp(1);
%             lowTime = bE(peakIdx)-0.001;
%             highTime = bE(peakIdx)+0.001;
        end
        baseline = sum(spikeTimeData{c,w} < min(0,lowTime))/(min(0,lowTime)+preTime); % spikes/ms
        ROI = sum(spikeTimeData{c,w} > lowTime & spikeTimeData{c,w} < highTime)/(highTime-lowTime); % spikes/ms
        prob(c,w) = ((ROI - baseline)*(highTime-lowTime))/sum(cds.waveforms.waveSent == w & cds.waveforms.chanSent == chanList(c));
    end
end

probOut = prob(chans,waveformTypes);

end