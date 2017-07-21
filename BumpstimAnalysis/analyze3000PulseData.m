%% set file names 

folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Mihili_20170712\';
% folderpath = 'D:\Lab\Data\StimArtifact\Mihili\20170713_stimRecord\';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left PMd SN 6251-001460.cmp';
mapFileName = 'C:\Users\Joseph\Desktop\Mihili Left PMd SN 6251-001460.cmp';

pwd=cd;
cd(folderpath)
fileList = dir('*_processed.mat');

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
figDir = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Mihili_20170712\Summary Figures\';
figPrefix = 'Mihili_20170712_';
saveFigures = 0;

nn = 3;
%%
plotRasterStim(cds,nn,'makeFigure',1,'makeSubplots',0,'plotTitle',1,'waveformTypes',[1:1:numel(cds.waveforms.parameters)],...
    'chans',[1:1:numel(unique(cds.waveforms.chanSent))],'preTime',10/1000,'postTime',60/1000,'plotSpikeWaveforms',1,'timeAfterStimRawNoStim',20/1000,...
    'timeAfterStimRawArtifact',5/1000,'plotArtifacts',1,'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
    'maxArtifactsPerPlot',5,'plotFiltered',0,'stimsPerTrain',1,'stimElectrode',stimElectrode);

%% plot grid
plotArrayMap(cds,nn,mapFileName,'numRows',10,'numCols',10,...
    'stimElectrode',[20,52,56],'stimElectrodeColor',{'k','r','b'},'stimElectrodeLabel','string',...
    'recordingElectrode',cds.units(nn).chan,'recordingElectrodeColor','k','recordingElectrodeLabel','string',...
    'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix)
% 
%% interspike interval histogram

plotInterspikeIntervalHistogram(cds,nn,'xLim',[0,20],'binSize',0.2,'displayText',1,...
    'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix);

%% probability of eliciting a response
nn=124;
probOut = getProbabilityOfResponse(cds,nn,'peakPeriod','automatic','preTime',20/1000,'postTime',60/1000,...
    'chans',[1:1:numel(unique(cds.waveforms.chanSent))],'waveformTypes',[1:1:numel(cds.waveforms.parameters)])

%% double peak dependence -- if applicable
% nn=99;
nn = 80;
[probActual,probExpected] = getDoublePeakProbability(cds,nn,'peakPeriod1',[3,4]/1000,'peakPeriod2',[4.3,5.5]/1000,'preTime',20/1000,'postTime',60/1000,...
    'chans',[1:1:numel(unique(cds.waveforms.chanSent))],'waveformTypes',[1:1:numel(cds.waveforms.parameters)])

%% plot artifacts
nn = 2;

plotArtifactsStim(cds,nn,1,1,'rowSubplot',3,'colSubplot',3,...
    'maxArtifactsPerPlot',8,'plotArtifactsSeparated',0,'plotTitle',0,...
    'plotFiltered',1,'randomSample',1,'templateSubtract',1,'plotXRange',[1,60])

%% plot PSTH
nn=124;
saveFigures = 0;
plotPSTHStim(cds,nn,'binSize',0.2/1000,'makeFigure',1,'makeSubplots',0,'plotTitle',1,'waveformTypes',[1:1:numel(cds.waveforms.parameters)],...
            'chans',[1:1:numel(unique(cds.waveforms.chanSent))],'preTime',10/1000,'postTime',60/1000,'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
            'plotLine',1,'plotAllOnOneFigure',0,'lineColor',{'k','r','b','g','m'})

%% whole array analysis

plotPSTHStimWholeArray(cds,mapFileName,'preTime',10/1000,'postTime',30/1000,'binSize',0.0002,...
    'waveformTypes',4,'chans',1,'plotLine',1,'lineWidth',0.5,'plotStimOn',1,'plotStimChan',1,'plotProbabilityText',1);




%% make all and save all
figDir = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Mihili_20170713\Summary Figures\';
figPrefix = 'Mihili_20170713_';
saveFigures = 1;
stimElectrode = [57];
% stimElectrode = [13,42,57,70];
for nn =107:size(cds.units,2)
    if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255)

        plotRasterStim(cds,nn,'makeFigure',1,'makeSubplots',0,'plotTitle',1,'waveformTypes',[1:1:numel(cds.waveforms.parameters)],...
            'chans',[1:1:numel(unique(cds.waveforms.chanSent))],'preTime',10/1000,'postTime',30/1000,'plotSpikeWaveforms',1,'timeAfterStimRawNoStim',20/1000,...
            'timeAfterStimRawArtifact',5/1000,'plotArtifacts',1,'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
            'maxArtifactsPerPlot',5,'plotFiltered',[0,1],'stimsPerTrain',1,'stimElectrode',stimElectrode);
        
        plotArrayMap(cds,nn,mapFileName,'numRows',10,'numCols',10,...
            'stimElectrode',stimElectrode,'stimElectrodeColor',{'k','r','b','g'},'stimElectrodeLabel','string',...
            'recordingElectrode',cds.units(nn).chan,'recordingElectrodeColor','k','recordingElectrodeLabel','string',...
            'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix)

        plotInterspikeIntervalHistogram(cds,nn,'xLim',[0,20],'binSize',0.2,'displayText',1,...
            'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix);
        
        plotPSTHStim(cds,nn,'binSize',0.2/1000,'makeFigure',1,'makeSubplots',0,'plotTitle',1,'waveformTypes',[1:1:numel(cds.waveforms.parameters)],...
            'chans',[1:1:numel(unique(cds.waveforms.chanSent))],'preTime',10/1000,'postTime',30/1000,'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix)

        
        close all
   end
end

