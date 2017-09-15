%% set file names 

folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20170614_bumpstim\';
% folderpath = 'D:\Lab\Data\StimArtifact\Han\bumpstim\20170614\';
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
neuronNumber = 11;
for chan = 1:numel(unique(cds.waveforms.chanSent))
    for wave = 1:numel(unique(cds.waveforms.waveSent))
        plotWaveformsStim(cds,neuronNumber,chan,wave,'timeAfterStimRawNoStim',5/1000,'timeAfterStimRawArtifact',5/1000,'plotFiltered',[0,1]);
    end
end
%% 2. compare non-stim and stim trials to see response long-term (raster -- will need to be able to specify time window and center)
neuronNumber = 116;
optsTask = [];
optsTask.TRIAL_LIST = {'ctrHoldBump','ctrHoldBumpStim'}; % set to trials used
optsTask.ZERO_MARKER = 'bumpTime';
optsTask.COMBINE = {'stimCode'};
optsTask.IGNORE = {'tgtDir'};
optsTask.PLOT_STIM_TIME = 1;

optsPlot.MARKER_STYLE = 'line';
optsPlot.Y_LABEL = 'Trial';
optsPlot.X_LABEL = 'Time after bump onset (s)';
optsPlot.X_LIMITS = [-0.5,1];
optsSave = [];

plotRasterBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);

optsTask.COMBINE = {'bumpDir'};
plotRasterBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);

% 
optsTask.TRIAL_LIST = {'noBump'};
optsTask.COMBINE = {'tgtDir'};
optsTask.IGNORE = {'bumpDir'};
optsTask.ZERO_MARKER = 'goCueTime';
optsPlot.Y_LABEL = 'Trial number';
optsPlot.X_LABEL = 'Time after go cue (s)';
plotRasterBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);

% across entire trial

%% 3. same as above as a PSTH/line PSTH. 

% neuronNumber = 11;
optsTask.TRIAL_LIST = {'ctrHoldBump','ctrHoldBumpStim'}; % set to trials used
optsTask.ZERO_MARKER = 'bumpTime';
optsTask.COMBINE = {'stimCode'};
optsTask.IGNORE = {'tgtDir'};
optsTask.BIN_SIZE = 20/1000;
optsTask.SAME_X_LIMITS = 1;
optsTask.SAME_Y_LIMITS = 1;

optsPlot.BAR_STYLE = 'line';
optsPlot.Y_LABEL = 'Average firing rate (spikes/s)';
optsPlot.X_LABEL = 'Time after bump onset (s)';
optsPlot.LEGEND_STRING = {'No stim','Stim'};
optsSave = [];
plotPSTHBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);

optsTask.COMBINE = {'bumpDir'};
optsPlot.LEGEND_STRING = {'Right','Down','Left','Up'};
plotPSTHBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);

%% raster for whole array
neuronNumber = 11;
optsArray.MAP_FILE = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

optsTask.TRIAL_LIST = {'ctrHoldBump'}; % set to trials used
optsTask.ZERO_MARKER = 'bumpTime';
optsTask.COMBINE = {'bumpDir'};
optsTask.IGNORE = {'tgtDir'};

optsPlot.MAKE_FIGURE = 0;

plotWholeArrayLIB(cds,'plotPSTHBumpstim',optsTask,optsPlot,optsSave,optsArray);

%% PSTH for whole array



% 4. assess ability to record directly after stimluation


%% plot all units and save to a folder

