%% set file names 

folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Mihili_20170721\';
% folderpath = 'D:\Lab\Data\StimArtifact\Mihili\20170713_stimRecord\';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left PMd SN 6251-001460.cmp';
% mapFileName = 'C:\Users\Joseph\Desktop\Mihili Left PMd SN 6251-001460.cmp';

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
    stimElectrode = 57;
end
load(fileList(fileNumber).name);
cd(pwd);

figDir = '';
figPrefix = '';
saveFigures = 0;

%% plot raster, waves, and PSTH for a give neuron number
figDir = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Mihili_20170712\Summary Figures\';
% figDir = strcat(cd,'\');
figPrefix = 'Mihili_20170712_';
saveFigures = 0;


%%
nn = 85;
plotWaveforms = 0;
plotArtifacts = 0;
plotRasterStim(cds,nn,'plotSpikeWaveforms',plotWaveforms,'plotArtifacts',plotArtifacts,'stimsPerTrain',1,...
    'waveformTypes',[1:1:numel(cds.waveforms.parameters)],'chans',[1:1:numel(unique(cds.waveforms.chanSent))],...
    'preTime',10/1000,'postTime',30/1000,'timeAfterStimRawNoStim',20/1000,'timeAfterStimRawArtifact',5/1000,...
    'markerstyle','dot','linelength',1,'linewidth',1.5,...
    'makeFigure',1,'makeSubplots',0,'plotTitle',1,...        
    'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
    'maxArtifactsPerPlot',5,'plotArtifactsFiltered',[0,1],'stimElectrode',stimElectrode, 'plotWaveformsFiltered',[0,1],...
    'maxWaveformsPlot',100,'plotStimuliGroup','all'); % plotStimuliGroup = 'responsive','nonResponsive','all'

%% plot grid
plotArrayMap(cds,nn,mapFileName,'numRows',10,'numCols',10,...
    'stimElectrode',[15,23,36,40,45],'stimElectrodeColor',{'k','r','b','g','m'},'stimElectrodeLabel','string',...
    'recordingElectrode',cds.units(nn).chan,'recordingElectrodeColor','k','recordingElectrodeLabel','string',...
    'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix)

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
nn = 37;

plotArtifactsStim(cds,nn,1,1,'rowSubplot',2,'colSubplot',1,...
    'maxArtifactsPerPlot',4,'plotArtifactsSeparated',1,'plotTitle',0,...
    'plotFiltered',[0],'randomSample',0,'templateSubtract',1,'plotXRange',[1,300])

plotArtifactsStim(cds,nn,1,1,'rowSubplot',2,'colSubplot',1,...
    'maxArtifactsPerPlot',4,'plotArtifactsSeparated',1,'plotTitle',0,...
    'plotFiltered',[1],'randomSample',0,'templateSubtract',1,'plotXRange',[1,300])

%% plot PSTH
nn=112;
saveFigures = 0;
groups = {'all','responsive','nonResponsive'};
plotPSTHStim(cds,nn,'binSize',0.2/1000,'makeFigure',1,'makeSubplots',0,'plotTitle',1,'waveformTypes',[1:1:numel(cds.waveforms.parameters)],...
        'chans',[1:1:numel(unique(cds.waveforms.chanSent))],'preTime',10/1000,'postTime',90/1000,'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
        'plotLine',0,'plotAllOnOneFigure',0,'lineColor',{'k','r','b','g','m'},...
        'plotStimuliGroup','responsive')

%% whole array analysis

plotPSTHStimWholeArray(cds,mapFileName,'preTime',10/1000,'postTime',30/1000,'binSize',0.0002,...
    'waveformTypes',4,'chans',1,'plotLine',1,'lineWidth',0.5,'plotStimOn',1,'plotStimChan',1,'plotProbabilityText',1);




%% make all and save all
figDir = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Mihili_20170720\Summary Figures\';
figPrefix = 'Mihili_20170720_';
saveFigures = 1;
% stimElectrode = [57];
stimElectrode = [31,44,49,52,96];
for nn =1:size(cds.units,2)
    if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255)

        plotRasterStim(cds,nn,'makeFigure',1,'makeSubplots',0,'plotTitle',1,'waveformTypes',[1:1:numel(cds.waveforms.parameters)],'markersize',4,...
            'chans',[1:1:numel(unique(cds.waveforms.chanSent))],'preTime',10/1000,'postTime',30/1000,'plotSpikeWaveforms',1,'timeAfterStimRawNoStim',20/1000,...
            'timeAfterStimRawArtifact',5/1000,'plotArtifacts',1,'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
            'maxArtifactsPerPlot',5,'plotFiltered',[0,1],'stimsPerTrain',1,'stimElectrode',stimElectrode,'plotWaveformsFiltered',[0,1]);
        
        plotArrayMap(cds,nn,mapFileName,'numRows',10,'numCols',10,...
            'stimElectrode',stimElectrode,'stimElectrodeColor',{'k','r','b','g','m'},'stimElectrodeLabel','string',...
            'recordingElectrode',cds.units(nn).chan,'recordingElectrodeColor','k','recordingElectrodeLabel','string',...
            'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix)

        plotInterspikeIntervalHistogram(cds,nn,'xLim',[0,20],'binSize',0.2,'displayText',1,...
            'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix);
        
        plotPSTHStim(cds,nn,'binSize',0.2/1000,'makeFigure',1,'makeSubplots',0,'plotTitle',1,'waveformTypes',[1:1:numel(cds.waveforms.parameters)],...
            'chans',[1:1:numel(unique(cds.waveforms.chanSent))],'preTime',10/1000,'postTime',30/1000,'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix)

        plotPSTHStim(cds,nn,'binSize',0.2/1000,'makeFigure',1,'makeSubplots',0,'plotTitle',0,'waveformTypes',[1:1:numel(cds.waveforms.parameters)],...
            'chans',[1:1:numel(unique(cds.waveforms.chanSent))],'preTime',10/1000,'postTime',30/1000,'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
            'plotLine',1,'plotAllOnOneFigure',1,'lineColor',{'k','r','b','g','m'})
        
        close all
   end
end

