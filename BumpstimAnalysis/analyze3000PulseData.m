%% set file names 
folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Chips_20171025\';
% folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Chips_20171024\';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left PMd SN 6251-001460.cmp';
mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
% mapFileName = 'C:\Users\Joseph\Desktop\Mihili Left PMd SN 6251-001460.cmp';

pwd=cd;
cd(folderpath)
fileList = dir('*all_processed*');

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

%% extract relevant data for a given unit
tic
optsExtract.NEURON_NUMBER = 3;

optsExtract.STIMULI_RESPONSE = 'all';
optsExtract.STIM_ELECTRODE = 50;
optsExtract.STIMULATIONS_PER_TRAIN = 1;
optsExtract.STIMULATION_BATCH_SIZE = 4000;

optsExtract.PRE_TIME = 10/1000; % made negative in the function
optsExtract.POST_TIME = 90/1000;
optsExtract.BIN_SIZE = 0.2/1000;
optsExtract.TIME_AFTER_STIMULATION_WAVEFORMS = 10/1000;

unitData = extractDataAroundStimulations(cds,optsExtract);
toc
%% plot raster, and PSTH for the given unit above
optsPlotFunc.BIN_SIZE = optsExtract.BIN_SIZE;

optsPlotFunc.PRE_TIME = 5/1000;
optsPlotFunc.POST_TIME = 30/1000;

rasterPlots = plotRasterStim(unitData,optsExtract.NEURON_NUMBER,optsPlotFunc);

optsPlotFunc.PLOT_ALL_ONE_FIGURE = 0;
optsPlotFunc.PLOT_LINE = 1;

optsPlotFunc.SMOOTH = 0;
optsPlotFunc.SMOOTH_STD_DEV = 0.35;

PSTHPlots = plotPSTHStim(unitData,optsExtract.NEURON_NUMBER,optsPlotFunc);

%% plot waveforms (raw and filtered)

optsWaveforms = [];

WaveformPlots = plotWaveformsStim(cds,optsExtract.NEURON_NUMBER,optsWaveforms);

%% plot raster, waves, and PSTH for a give neuron number
figDir = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Chips_20171025\';
% figDir = strcat(cd,'\');
figPrefix = 'Chips_20171025_';
saveFigures = 0;


%%
nn = 81;
plotWaveforms = 0;
plotArtifacts = 1;
outputData = plotRasterStim(cds,nn,'plotSpikeWaveforms',plotWaveforms,'plotArtifacts',plotArtifacts,'stimsPerTrain',1,...
    'waveformTypes',[1:1:numel(unique(cds.waveforms.waveSent))],'chans',[1:1:numel(unique(cds.waveforms.chanSent))],...
    'preTime',10/1000,'postTime',110/1000,'timeAfterStimRawNoStim',10/1000,'timeAfterStimRawArtifact',5/1000,...
    'markerstyle','.','markersize',4,'linelength',1,'linewidth',1.5,...
    'makeFigure',1,'makeSubplots',0,'plotTitle',0,'plotXRange',[0,5],...        
    'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
    'maxArtifactsPerPlot',30,'plotArtifactsFiltered',[1],'stimElectrode',stimElectrode, 'plotWaveformsFiltered',[0,1],...
    'maxWaveformsPlot',100,'plotStimuliGroup','all','sortData','no'); % plotStimuliGroup = 'responsive','nonResponsive','all'

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
% nn = 112;
nn=67;
figDir = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Chips_20171024_dukeBoardchan65\';
% figDir = strcat(cd,'\');
figPrefix = 'Chips_20171024_dukeBoardV1_';
saveFigures = 0;

for wave = 1:numel(unique(cds.waveforms.waveSent))  
    plotArtifactsStim(cds,nn,1,wave,'rowSubplot',1,'colSubplot',1,...
        'maxArtifactsPerPlot',10,'plotArtifactsSeparated',0,'plotTitle',0,...
        'plotFiltered',[0],'randomSample',1,'templateSubtract',0,'plotYRange',[-9000,9000],'plotXRange',[0,8],...
        'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
        'divisor',1)
    
%     plotArtifactsStim(cds,nn,1,wave,'rowSubplot',1,'colSubplot',1,...
%         'maxArtifactsPerPlot',10,'plotArtifactsSeparated',0,'plotTitle',0,...
%         'plotFiltered',[1],'randomSample',1,'templateSubtract',0,'plotYRange',[-400,400],'plotXRange',[0,8],...
%         'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
%         'divisor',1)
end
% plotArtifactsStim(cds,nn,find(unique(cds.waveforms.chanSent)==45),1,'rowSubplot',1,'colSubplot',1,...
%     'maxArtifactsPerPlot',4,'plotArtifactsSeparated',1,'plotTitle',0,...
%     'chans',[find(unique(cds.waveforms.chanSent)==45)],...
%     'plotFiltered',[1],'randomSample',0,'templateSubtract',0,'plotXRange',[1,300])

%% plot PSTH
nn=12;
saveFigures = 0;
groups = {'all','responsive','nonResponsive'};
plotPSTHStim(cds,nn,'binSize',0.2/1000,'makeFigure',1,'makeSubplots',0,'plotTitle',1,'stimsPerTrain',1,'stimCenter',1,...
    'waveformTypes',[1:1:numel(unique(cds.waveforms.waveSent))],'chans',[1:1:numel(unique(cds.waveforms.chanSent))],...
    'preTime',10/1000,'postTime',30/1000,'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
    'plotLine',1,'plotAllOnOneFigure',1,'lineColor',{'k','r','b',[0 0.5 0],'m',[0.5,0.5,0.5],[1,0.6,0]},...
    'plotStimuliGroup','nonresponsive');

%% whole array plot

plotPSTHStimWholeArray(cds,mapFileName,'preTime',10/1000,'postTime',30/1000,'binSize',0.0002,...
    'waveformTypes',[1:numel(unique(cds.waveforms.waveSent))],'chans',[1:numel(unique(cds.waveforms.chanSent))],'plotLine',1,'lineWidth',0.5,'plotStimOn',1,'plotStimChan',1,'plotProbabilityText',1);

%% heatmap
% for i = 0:0.5:5.5
    plotHeatmap(cds,mapFileName,'baselinePreTime',-10/1000,'baselinePostTime',-2/1000,...
        'stimPreTime',1.5/1000,'stimPostTime',5/1000,...
        'binSize',0.0002,'waveformTypes',1:numel(unique(cds.waveforms.waveSent)),'chans',1:numel(unique(cds.waveforms.chanSent)),'plotLine',1,...
        'lineWidth',0.5,'plotStimOn',1,'plotStimChan',1,...
        'saveFigures',saveFigures,'figDir',figDir,'figPrefix',strcat(figPrefix));
    
%     close all
% end
%% make all and save all

figDir = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Chips_20171026\';
figPrefix = 'Chips_20171026_';

saveFigures = 1;
% stimElectrode = [57];
stimElectrode = unique(cds.waveforms.chanSent);
plotWaveforms = 1;
plotArtifacts = 0;

plotHeatmap(cds,mapFileName,'baselinePreTime',-10/1000,'baselinePostTime',-2/1000,...
    'stimPreTime',0,'stimPostTime',5/1000,...
    'binSize',0.0002,'waveformTypes',1:numel(unique(cds.waveforms.waveSent)),'chans',1:numel(unique(cds.waveforms.chanSent)),'plotLine',1,...
    'lineWidth',0.5,'plotStimOn',1,'plotStimChan',1,'plotProbabilityText',1,...
    'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix);

% for nn = 1:size(cds.units,2)
%     if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255)
% % 
%         plotRasterStim(cds,nn,'plotSpikeWaveforms',plotWaveforms,'plotArtifacts',plotArtifacts,'stimsPerTrain',1,...
%             'waveformTypes',[1:1:numel(unique(cds.waveforms.waveSent))],'chans',[1:1:numel(unique(cds.waveforms.chanSent))],...
%             'preTime',10/1000,'postTime',30/1000,'timeAfterStimRawNoStim',10/1000,'timeAfterStimRawArtifact',8/1000,...
%             'markerstyle','.','linelength',1,'linewidth',1.5,...
%             'makeFigure',1,'makeSubplots',0,'plotTitle',0,...        
%             'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
%             'maxArtifactsPerPlot',10,'plotArtifactsFiltered',[0,1],'stimElectrode',stimElectrode, 'plotWaveformsFiltered',[0,1],...
%             'maxWaveformsPlot',100,'plotStimuliGroup','all','sortData','no'); % plotStimuliGroup = 'responsive','nonResponsive','all'
% 
%         plotArrayMap(cds,nn,mapFileName,'numRows',10,'numCols',10,...
%             'stimElectrode',stimElectrode,'stimElectrodeColor',{[0 0.5 0],'r','b','k','m',[0.5,0.5,0.5],[1,0.6,0]},'stimElectrodeLabel','string',...
%             'recordingElectrode',cds.units(nn).chan,'recordingElectrodeColor','k','recordingElectrodeLabel','string',...
%             'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix)
% 
%         plotInterspikeIntervalHistogram(cds,nn,'xLim',[0,20],'binSize',0.2,'displayText',1,...
%             'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix);
%         
%         plotPSTHStim(cds,nn,'binSize',0.2/1000,'makeFigure',1,'makeSubplots',0,'plotTitle',0,'waveformTypes',[1:1:numel(cds.waveforms.parameters)],...
%             'chans',[1:1:numel(unique(cds.waveforms.chanSent))],'preTime',10/1000,'postTime',30/1000,'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
%             'barColor','k')
% 
%         plotPSTHStim(cds,nn,'binSize',0.2/1000,'makeFigure',1,'makeSubplots',0,'plotTitle',0,'waveformTypes',[1:1:numel(cds.waveforms.parameters)],...
%             'chans',[1:1:numel(unique(cds.waveforms.chanSent))],'preTime',10/1000,'postTime',30/1000,'saveFigures',saveFigures,'figDir',figDir,'figPrefix',figPrefix,...
%             'plotLine',1,'plotAllOnOneFigure',1,'lineColor',{[0 0.5 0],'r','b','k','m',[0.5,0.5,0.5],[1,0.6,0]})
%         
%         close all
%    end
% end

