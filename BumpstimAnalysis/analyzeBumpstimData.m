%% set file names 

folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20170614_bumpstim\';
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
    stimElectrode = 42;
end
load(fileList(fileNumber).name);
cd(pwd);

figDir = '';
figPrefix = '';
saveFigures = 0;

%% 1. look waveforms around stimulation and not around stimulation to
% confirm that it is a neuron
neuronNumber = 93;
for chan = 1:numel(unique(cds.waveforms.chanSent))
    for wave = 1:numel(unique(cds.waveforms.waveSent))
        plotWaveformsStim(cds,neuronNumber,chan,wave,'timeAfterStimRawNoStim',5/1000,'timeAfterStimRawArtifact',5/1000,'plotFiltered',[0,1]);
    end
end
%% 2. compare non-stim and stim trials to see response long-term (raster -- will need to be able to specify time window and center)
neuronNumber = 93;
optsTask = [];
% optsTask.TRIAL_LIST = {'ctrHoldBump';'ctrHoldBumpStim'}; % set to trials used
optsPlot = [];
optsSave = [];

plotRasterBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);
% across entire trial
%% 3. same as above as a PSTH/line PSTH. 

neuronNumber = 101;
optsTask.TRIAL_LIST = {'ctrHoldBump';'ctrHoldBumpStim'}; % set to trials used
optsPlot = [];
optsSave = [];

plotRasterBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);
% 4. assess ability to record directly after stimluation

