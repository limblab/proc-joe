%% set file names 
folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20170628\';
pwd=cd;
cd(folderpath)
fileList = dir('*_processed.mat');

% load file
load(fileList(1).name);
cd(pwd);

%% plot raster, waves, and PSTH for a give neuron number
figDir = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20170628\Summary Figures\';
figPrefix = 'Han_20170628_chan56stim_30uA_';
saveFigures = 0;

nn = 101;

plotRaster(cds,nn);

%% Raster for a given index in cds.units -- basic analysis

    
nn = 101;

timeBeforeStim = 10/1000; % 10 ms in seconds
timeAfterStim = 20/1000; % 10 ms in seconds
figure;
for st = 1:numel(cds.stimOn)-1
    spikeMask = cds.units(nn).spikes.ts > cds.stimOn(st)-timeBeforeStim & cds.units(nn).spikes.ts < cds.stimOn(st)+timeAfterStim;
    spikesPlot = (cds.units(nn).spikes.ts(spikeMask) - cds.stimOn(st))*1000;
    if(~isempty(spikesPlot))
        plot(spikesPlot,st,'k.')
    end
    xlim([-1*timeBeforeStim*1000 timeAfterStim*1000])
    hold on
end    

ylabel('Stimuli')
xlabel('Time after stimulation (ms)')
ylim([0 3000])
formatForLee(gcf);
figName = strcat(figPrefix,'_chan',num2str(cds.units(nn).chan),'_raster');
saveFigure(gcf,figDir,figName);
% plot all neuron waves
f=figure
maxWavesPlot = 100;

subplot(3,1,1) % non stim region waves
timeAfterStim = 20/1000;
wavesPlot = [];
spikeMask = zeros(numel(cds.units(nn).spikes.ts,1));
for st = 1:numel(cds.stimOn)
    if(st==numel(cds.stimOn))
        spikeMask = spikeMask | (cds.units(nn).spikes.ts > cds.stimOn(st) + timeAfterStim);
    else
        spikeMask = spikeMask | (cds.units(nn).spikes.ts > cds.stimOn(st) + timeAfterStim & cds.units(nn).spikes.ts < cds.stimOn(st+1));
    end
end

waveIdx = find(spikeMask);
waveIdx = datasample(waveIdx,min(maxWavesPlot,numel(waveIdx)),'Replace',false);
if(numel(waveIdx) > 0)
    for wave = 1:numel(waveIdx)
        rawIdx = getRawDataIdx(cds.units(nn).spikes.ts(waveIdx(wave)),cds.units(nn).chan,cds.rawData.ts,cds.rawData.elec);
        if(rawIdx ~= -1)
            wavesPlot = [wavesPlot;cds.rawData.waveforms(rawIdx,:)];
        end
    end
    xDataWaves = ((1:size(wavesPlot,2))-1)/30; % in ms
    plot(xDataWaves,wavesPlot)
end

ylim([-300 300])

ylabel('Voltage (\muV)')
xlabel('Time (ms)')
formatForLee(gcf);


subplot(3,1,2) % stimulation artifacts aligned to waveform
timeAfterStim = 4/1000;
timeBeforeStim = 0/1000;
wavesPlot = [];
spikeMask = zeros(numel(cds.units(nn).spikes.ts,1));
artifactMask = zeros(numel(cds.stimOn),1);
for st = 1:numel(cds.stimOn)
    spikeMask = spikeMask | (cds.units(nn).spikes.ts < cds.stimOn(st) + timeAfterStim & cds.units(nn).spikes.ts > cds.stimOn(st) + timeBeforeStim);
    if(sum((cds.units(nn).spikes.ts < cds.stimOn(st) + timeAfterStim & cds.units(nn).spikes.ts > cds.stimOn(st)))>0)
        artifactMask(st,1) = 1;
    end
end

waveIdx = find(spikeMask);
waveIdx = datasample(waveIdx,min(maxWavesPlot,numel(waveIdx)),'Replace',false);
if(numel(waveIdx) > 0)
    for wave = 1:numel(waveIdx)
        rawIdx = getRawDataIdx(cds.units(nn).spikes.ts(waveIdx(wave)),cds.units(nn).chan,cds.rawData.ts,cds.rawData.elec);
        if(rawIdx ~= -1)
            wavesPlot = [wavesPlot;cds.rawData.waveforms(rawIdx,:)];
        end
    end
    xDataWaves = ((1:size(wavesPlot,2))-1)/30; % in ms
    plot(xDataWaves,wavesPlot)
end

ylim([-300 300])
ylabel('Voltage (\muV)')
xlabel('Time (ms)')
formatForLee(gcf);


subplot(3,1,3) % aligned to stimulation artifact
timeAfterArtifact = 8/1000;
chan = cds.units(nn).chan;
artifactsPlot = squeeze(cds.artifactData.artifact(artifactMask==1,chan,1:timeAfterArtifact*30000));
xDataArtifact = ((1:size(artifactsPlot,2))-1)/30;

plot(xDataArtifact,artifactsPlot)
ylim([-400 400])
ylabel('Voltage (\muV)')
xlabel('Time (ms)')
formatForLee(gcf);

figName = strcat(figPrefix,'_chan',num2str(cds.units(nn).chan),'_waves');
saveFigure(f,figDir,figName);

% PSTH
timeBeforeStim = 10/1000; % 10 ms in seconds
timeAfterStim = 40/1000; % 10 ms in seconds
binSize = 0.1/1000; % in seconds
figure;
for st = 1:numel(cds.stimOn)-1
    spikeMask = cds.units(nn).spikes.ts > cds.stimOn(st)-timeBeforeStim & cds.units(nn).spikes.ts < cds.stimOn(st)+timeAfterStim;
    spikesPlot = [spikesPlot;(cds.units(nn).spikes.ts(spikeMask) - cds.stimOn(st))*1000];

    xlim([-1*timeBeforeStim*1000 timeAfterStim*1000])
    hold on
end  

bE = (-1*timeBeforeStim:binSize:timeAfterStim)*1000;
[bC,bE] = histcounts(spikesPlot,bE);
bC = bC;
bar(bE(1:end-1) + (bE(2)-bE(1))/2, bC);

ylabel('Spike count')
xlabel('Time after stimulation (ms)')
formatForLee(gcf);

figName = strcat(figPrefix,'_chan',num2str(cds.units(nn).chan),'_PSTH');
saveFigure(gcf,figDir,figName);

% latency bar plot (zoomed in more PSTH basically)
figure;
timeBeforeStim = -1/1000; % 10 ms in seconds
timeAfterStim = 4/1000; % 10 ms in seconds
binSize = 0.1/1000; % in seconds
bE = (-1*timeBeforeStim:binSize:timeAfterStim)*1000;
[bC,bE] = histcounts(spikesPlot,bE);
bC = bC;
bar(bE(1:end-1) + (bE(2)-bE(1))/2, bC);

ylabel('Spike count')
xlabel('Time after stimulation (ms)')
formatForLee(gcf);

figName = strcat(figPrefix,'_chan',num2str(cds.units(nn).chan),'_PSTHzoomed');
saveFigure(gcf,figDir,figName);

% latency vs time since last spike
figure;
timeAfterStim = 3/1000;
wavesPlot = [];
spikeMask = zeros(numel(cds.units(nn).spikes.ts,1));
latency = [];
timeBetween = [];
for st = 1:numel(cds.stimOn)
    spikeMask = (cds.units(nn).spikes.ts < cds.stimOn(st) + timeAfterStim & cds.units(nn).spikes.ts > cds.stimOn(st));
    spikeIdx = find(spikeMask);
    spikeTs = cds.units(nn).spikes.ts(spikeIdx) - cds.stimOn(st);
    prevSpikeTs = cds.units(nn).spikes.ts(spikeIdx-1)-cds.stimOn(st);
    latency(end+1:end+numel(spikeTs)) = spikeTs*1000;
    timeBetween(end+1:end+numel(prevSpikeTs)) = (spikeTs - prevSpikeTs)*1000;
end

plot(timeBetween,latency,'k.')
ylabel('Latency (ms)')
xlabel('Time since last action potential (ms)')

figName = strcat(figPrefix,'_chan',num2str(cds.units(nn).chan),'_PostStimInterval');
saveFigure(gcf,figDir,figName);

%% For a unit -- probability of eliciting a spike
nn = 9;
timeAfterStim = 5/1000;
numEvoked = 0;
for st = 1:numel(cds.stimOn)
    spikeMask = (cds.units(nn).spikes.ts < cds.stimOn(st) + timeAfterStim & cds.units(nn).spikes.ts > cds.stimOn(st));
    if(sum(spikeMask) > 0)
        numEvoked = numEvoked + 1;
    end
end


