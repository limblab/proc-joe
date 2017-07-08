%% set file names 
folderpath = 'D:\Lab\Data\StimArtifact\Han\20170706_stimRecord\';
mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left PMd SN 6251-001460.cmp';
pwd=cd;
cd(folderpath)
fileList = dir('*all_processed.mat');

%% load file and parse for stim electrode number
fileNumber = 1;
chanIdx = strfind(fileList(fileNumber).name,'chan');
stimIdx = strfind(fileList(fileNumber).name,'stim');
if(~isempty(chanIdx) && numel(chanIdx) == 1 && ~isempty(stimIdx) && numel(stimIdx) == 1)
    stimElectrode = str2num(fileList(fileNumber).name(chanIdx+4:stimIdx-1));
else % manually input stim electrode
    stimElectrode = 17;
end
load(fileList(fileNumber).name);
cd(pwd);

%% plot raster, waves, and PSTH for a give neuron number
figDir = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20170629\Summary Figures\';
figPrefix = 'Han_20170628_chan42stim_250us';
saveFigures = 0;

nn = 108;

plotRasterStim(cds,nn,'makeFigure',1,'makeSubplots',0,'plotTitle',1,'waveformTypes',[1:1:numel(cds.waveforms.parameters)],...
    'preTime',10/1000,'postTime',30/1000,'plotSpikeWaveforms',1,'timeAfterStimRawNoStim',20/1000,...
    'timeAfterStimRawArtifact',9/1000,'plotArtifacts',1,'saveFigure',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
    'maxArtifactsPerPlot',5,'plotFiltered',0);

% plot grid
plotArrayMap(cds,nn,mapFileName,'numRows',10,'numCols',10,...
    'stimElectrode',stimElectrode,'stimElectrodeColor','k','stimElectrodeLabel','string',...
    'recordingElectrode',cds.units(nn).chan,'recordingElectrodeColor','k','recordingElectrodeLabel','string')

plotInterspikeIntervalHistogram(cds,nn,'xLim',[0,20],'binSize',0.2,'displayText',1);

%% plot PSTH
saveFigures = 0;

nn = 108;

plotPSTHStim(cds,nn,'binSize',0.2/1000,'makeFigure',1,'makeSubplots',0,'plotTitle',1,'waveformTypes',[1:1:numel(cds.waveforms.parameters)],...
    'preTime',10/1000,'postTime',20/1000,'saveFigure',saveFigures,'figDir',figDir,'figPrefix',figPrefix)

% plotLatencyVsSpikeTiming

% probability of eliciting a spike

% whole array analysis


%% Raster for a given index in cds.units -- basic analysis
figDir = 'D:\Lab\Data\StimArtifact\Han\20170622_3000pulses\Summary Figures\chan51stim_40uA_300us\';
figPrefix = 'Han_20170621_chan51stim_40uA_300us_';
    
nn = 9;

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
%%
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

for nn = 1:size(cds.units,2)
    if(cds.units(nn).ID ~= 255 && cds.units(nn).ID~=0)
        timeBeforeStim = 10/1000; % 10 ms in seconds
        timeAfterStim = 30/1000; % 10 ms in seconds
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
        if(size(artifactsPlot,1)>0)
            plot(xDataArtifact,artifactsPlot)
            ylim([-400 400])
            ylabel('Voltage (\muV)')
            xlabel('Time (ms)')
            formatForLee(gcf);

            figName = strcat(figPrefix,'_chan',num2str(cds.units(nn).chan),'_waves');
            saveFigure(f,figDir,figName);
        end

        % PSTH
        timeBeforeStim = 10/1000; % 10 ms in seconds
        timeAfterStim = 30/1000; % 10 ms in seconds
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
        timeAfterStim = 8/1000; % 10 ms in seconds
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
        close all
    end
end
%% probabilty of eliciting a spike
cd(folderpath)
mapFile = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
arrayMap = loadMapFile(mapFile);
timeAfterStim = 5/1000;
numEvoked = zeros(11,4);
arrayPos = zeros(11,2);

for f = 1:5
    counter = 1;
    load(fileList(f).name);

    for nn = 1:size(cds.units,2)
        if(cds.units(nn).ID ~= 255 && cds.units(nn).ID~=0)
            temp = 0;
            baseline = 0;
            for st = 1:numel(cds.stimOn)
                spikeMask = (cds.units(nn).spikes.ts < cds.stimOn(st) + timeAfterStim & cds.units(nn).spikes.ts > cds.stimOn(st));
                if(sum(spikeMask)>0)
                    temp = temp+1;
                end
            end
            chan(counter,1) = cds.units(nn).chan;
            numEvoked(counter,f) = temp/numel(cds.stimOn);
            arrayMapIdx = find(strcmp(cds.units(nn).label,arrayMap.label));
            arrayPos(counter,:) = [arrayMap.row(arrayMapIdx),arrayMap.col(arrayMapIdx)];
            counter=counter+1;
        end
    end

end


