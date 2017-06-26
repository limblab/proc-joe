% %% set file names 
% folderpath = 'R:\data\Han_13B1\Raw\20170622_3000pulseTrials\';
% pwd=cd;
% cd(folderpath)
% fileList = dir('*_cds_processed.mat');
% 
% %% load a cds
% load(fileList(1).name);
% cd(pwd);

%% Raster for a given index in cds.units
nn = 11;
    
timeBeforeStim = 10/1000; % 10 ms in seconds
timeAfterStim = 10/1000; % 10 ms in seconds
figure;
for st = 1:numel(cds.stimOn)
    spikeMask = cds.units(nn).spikes.ts-duration > cds.stimOn(st)-timeBeforeStim & cds.units(nn).spikes.ts-duration < cds.stimOn(st)+timeAfterStim;
    spikesPlot = (cds.units(nn).spikes.ts(spikeMask) - cds.stimOn(st) - duration)*1000;
    if(~isempty(spikesPlot))
        plot(spikesPlot,st,'k.')
    end
    xlim([-10 10])
    hold on
end    

% plot all neuron waves
figure
maxWavesPlot = 100;

subplot(3,1,1) % non stim region waves
timeAfterStim = 20/1000;
wavesPlot = [];
for st = 1:numel(cds.stimOn)
    if(st==numel(cds.stimOn))
        spikeMask = spikeMask | (cds.units(nn).spikes.ts-duration > cds.stimOn(st) + timeAfterStim);
    else
        spikeMask = spikeMask | (cds.units(nn).spikes.ts-duration > cds.stimOn(st) + timeAfterStim & cds.units(nn).spikes.ts-duration < cds.stimOn(st+1));
    end
end

waveIdx = find(spikeMask);
waveIdx = datasample(waveIdx,min(maxWavesPlot,numel(waveIdx)),'Replace',false);
if(numel(waveIdx) > 0)
    for wave = 1:numel(waveIdx)
        rawIdx = getRawDataIdx(cds.units(nn).spikes.ts(waveIdx(wave)),cds.units(nn).chan,cds.rawData.ts,cds.rawData.elec);
        if(rawIdx ~= -1)
            wavesPlot = [wavesPlot;rawData.waveforms(rawIdx,:)];
        end
    end
    xDataWaves = ((1:size(wavesPlot,2))-1)/30; % in ms
    plot(xDataWaves,wavesPlot(waveIdx,:))
end

% subplot(3,1,2) % stimulation artifacts aligned to waveform
% timeAfterStim = 2/1000;
% wavesPlot = [];
% for st = 1:numel(cds.stimOn)
%     spikeMask = cds.units(nn).spikes.ts-duration > cds.stimOn(st) & cds.units(nn).spikes.ts-duration < cds.stimOn(st) + timeAfterStim;
%     wavesPlot = [wavesPlot;cds.units(nn).spikes{spikeMask,2:end}];
% end
% 
% if(numel(waveIdx) > 0)
%     xDataWaves = ((1:size(wavesPlot,2))-1)/30; % in ms
%     plot(xDataWaves,wavesPlot(:,:))
% end
% 
% 
% subplot(3,1,3) % aligned to stimulation artifact

