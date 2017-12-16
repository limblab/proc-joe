%% set file names 

folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Chips_20171116_COstim\';
% folderpath = 'D:\Lab\Data\StimArtifact\Han\bumpstim\20170614\';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left PMd SN 6251-001460.cmp';
mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';

pwd=cd;
cd(folderpath)
fileList = dir('*_processed.mat');

%% load file and parse for stim electrode number
fileNumber = 2;
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
neuronNumber = 37;
for chan = 1:numel(unique(cds.waveforms.chanSent))
    for wave = 1:numel(unique(cds.waveforms.waveSent))
        plotWaveformsStim(cds,neuronNumber,chan,wave,'timeAfterStimRawNoStim',5/1000,'timeAfterStimRawArtifact',5/1000,'plotFiltered',[0,1]);
    end
end
%% 2. compare non-stim and stim trials to see response long-term (raster -- will need to be able to specify time window and center)
% neuronNumber = 95;

optsSave.FIGURE_SAVE = 0;
optsSave.FIGURE_DIR = 'D:\Lab\Data\StimArtifact\Han\bumpstim\20170614\';
optsSave.FIGURE_NAME_PREFIX = 'Han_20170614_CObump_chan42stim_';

optsTask = [];
optsTask.TRIAL_LIST = {'noStim','ctrHoldBumpStim'}; % set to trials used
optsTask.ZERO_MARKER = 'goCueTime';
% optsTask.STIM_CODE = [-1];
optsTask.COMBINE = {'tgtDir'};
optsTask.IGNORE = {'stimCode'};
optsTask.PLOT_STIM_TIME = 1;

optsPlot.MARKER_STYLE = 'line';
optsPlot.Y_LABEL = 'Trial';
optsPlot.X_LABEL = 'Time (s)';
optsPlot.X_LIMITS = [-1,2];

plotRasterBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);
%%
optsTask.COMBINE = {'bumpDir'};
plotRasterBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);
%%
optsTask.TRIAL_LIST = {'noBump'};
optsTask.COMBINE = {'tgtDir'};
optsTask.IGNORE = {'bumpDir'};
optsTask.ZERO_MARKER = 'goCueTime';
optsPlot.Y_LABEL = 'Trial number';
optsPlot.X_LABEL = 'Time after go cue (s)';
plotRasterBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);

%% 3. same as above as a PSTH/line PSTH. 

% neuronNumber = 11;
%%
neuronNumber = 116;
optsTask.TRIAL_LIST = {'ctrHoldBump','ctrHoldBumpStim'}; % set to trials used
optsTask.ZERO_MARKER = 'bumpTime';
optsTask.COMBINE = {'bumpDir'};
optsTask.COMBINE_SUBTRACT = 1;
optsTask.IGNORE = {'tgtDir'};
optsTask.BIN_SIZE = 20/1000;
optsTask.SAME_X_LIMITS = 1;
optsTask.SAME_Y_LIMITS = 1;

optsPlot.BAR_STYLE = 'line';
optsPlot.Y_LABEL = 'Average firing rate (spikes/s)';
optsPlot.X_LABEL = 'Time after bump onset (s)';
optsPlot.LEGEND_STRING = {'No stim','Stim'};

plotPSTHBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);
%%
optsTask.COMBINE = {'bumpDir'};
optsPlot.LEGEND_STRING = {'Right','Down','Left','Up'};
plotPSTHBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);

%% PSTH for whole array
optsArray.MAP_FILE = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
optsArray.STIM_ELEC = stimElectrode;

optsTask.TRIAL_LIST = {'ctrHoldBump','ctrHoldBumpStim'}; % set to trials used
optsTask.ZERO_MARKER = 'bumpTime';
optsTask.COMBINE = {'bumpDir'};
optsTask.COMBINE_SUBTRACT = 1;
optsTask.IGNORE = {'tgtDir'};

optsPlot.MAKE_FIGURE = 0;
optsPlot.LINE_WIDTH = 0.2;
plotWholeArrayLIB(cds,'plotPSTHBumpstim',optsTask,optsPlot,optsSave,optsArray);

% 4. assess ability to record directly after stimluation


%% plot all units and save to a folder

for neuronNumber = 35:size(cds.units,2)
    if(cds.units(neuronNumber).ID ~= 0 && cds.units(neuronNumber).ID ~= 255)
        try
            optsSave.FIGURE_SAVE = 1;
            optsSave.FIGURE_DIR = 'D:\Lab\Data\StimArtifact\Han\bumpstim\20170614\';
            optsSave.FIGURE_NAME_PREFIX = 'Han_20170614_CObump_chan42stim_';

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
            
            plotRasterBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);
            close all

            optsTask.COMBINE = {'bumpDir'};
            plotRasterBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);
            close all
            % 
            optsTask.TRIAL_LIST = {'noBump'};
            optsTask.COMBINE = {'tgtDir'};
            optsTask.IGNORE = {'bumpDir'};
            optsTask.ZERO_MARKER = 'goCueTime';
            optsPlot.Y_LABEL = 'Trial number';
            optsPlot.X_LABEL = 'Time after go cue (s)';
            plotRasterBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);
            close all
            % across entire trial

            % 3. same as above as a PSTH/line PSTH. 

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

            plotPSTHBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);
            close all
            optsTask.COMBINE = {'bumpDir'};
            optsPlot.LEGEND_STRING = {'Right','Down','Left','Up'};
            plotPSTHBumpstim(cds,neuronNumber,optsTask,optsPlot,optsSave);
            close all
        catch
        end
    end
end
