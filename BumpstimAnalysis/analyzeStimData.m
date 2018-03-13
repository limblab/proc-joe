%% set file names 
% folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Chips_20171025\';
folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Chips_20171024_dukeBoardchan65\';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left PMd SN 6251-001460.cmp';
mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';

pwd=cd;
cd(folderpath)
fileList = dir('*proc*');

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

optsExtract.NEURON_NUMBER = 81;

optsExtract.STIMULI_RESPONSE = 'all';
optsExtract.STIM_ELECTRODE = unique(cds.waveforms.chanSent);
optsExtract.STIMULATIONS_PER_TRAIN = 1;
optsExtract.STIMULATION_BATCH_SIZE = 1000;

optsExtract.PRE_TIME = 20/1000; % made negative in the function
optsExtract.POST_TIME = 180/1000;
optsExtract.BIN_SIZE = 0.2/1000;
optsExtract.TIME_AFTER_STIMULATION_WAVEFORMS = 10/1000;

unitData = extractDataAroundStimulations(cds,optsExtract);

%% plot raster, and PSTH for the given unit above
optsPlotFunc.BIN_SIZE = optsExtract.BIN_SIZE;
optsPlotFunc.FIGURE_SAVE = 1;
optsPlotFunc.FIGURE_DIR = folderpath;
optsPlotFunc.FIGURE_PREFIX = 'Chips_20171024_dukeBoardV1_velocityEffect';

optsPlotFunc.PRE_TIME = 20/1000;
optsPlotFunc.POST_TIME = 180/1000;
optsPlotFunc.SORT_DATA = 'velocity';

rasterPlots = plotRasterStim(unitData,optsExtract.NEURON_NUMBER,optsPlotFunc);

% optsPlotFunc.PLOT_ALL_ONE_FIGURE = 1;
% optsPlotFunc.PLOT_LINE = 1;
% optsPlotFunc.PLOT_TITLE = 0;
% 
% optsPlotFunc.SMOOTH = 0;
% optsPlotFunc.SMOOTH_STD_DEV = 0.35;
% optsPlotFunc.MAKE_LEGEND = 0;
% optsPlotFunc.LEGEND_STR = {'20\muA','40\muA'};
% PSTHPlots = plotPSTHStim(unitData,optsExtract.NEURON_NUMBER,optsPlotFunc);

%% plot waveforms (raw and filtered)
optsWaveforms = [];

optsWaveforms.FIGURE_SAVE = 0;
optsWaveforms.FIGURE_DIR = folderpath;
optsWaveforms.FIGURE_PREFIX = 'Chips_20171024_dukeBoardV1_';

optsWaveforms.YLIM = [-800,800];
optsWaveforms.TIME_AFTER_STIMULATION_ARTIFACT = 4/1000;
WaveformPlots = plotWaveformsStim(cds,optsExtract.NEURON_NUMBER,optsWaveforms);

%% plot artifacts (raw and filtered)

optsArtifacts = [];
optsArtifacts.WAVEFORMS_TO_PLOT = [4];
optsArtifacts.CHANS_TO_PLOT = [1];
optsArtifacts.MAX_WAVES_PLOT = 20;
optsArtifacts.YLIM = [-600,600];
optsArtifacts.XLIM = [0,10];
optsArtifacts.ADJUST_FOR_BIT_ERROR = 0;
optsArtifacts.ROW_SUBPLOT = 1;
optsArtifacts.COL_SUBPLOT = 1;
optsArtifacts.NUM_PAD = 200;
optsArtifacts.RANDOM_SAMPLE = 1;
optsArtifacts.ARTIFACT_MULTIPLIER = 5;

ArtifactPlots = plotArtifactsStim(cds,optsExtract.NEURON_NUMBER,optsArtifacts);

%% plot grid
optsGrid.STIM_ELECTRODE = unique(cds.waveforms.chanSent);
optsGrid.RECORDING_ELECTRODE = cds.units(optsExtract.NEURON_NUMBER).chan;

optsGrid.STIM_ELECTRODE_COLOR = {'k','r','b',[0,0.5,0],'m'};
optsGrid.STIM_ELECTRODE_LABEL = 'string';

optsGrid.RECORDING_ELECTRODE_COLOR = 'k';
optsGrid.RECORDING_ELECTRODE_LABEL = 'string';

ArrayPlots = plotArrayMap(cds,optsExtract.NEURON_NUMBER,mapFileName,optsGrid);

%%  ISI
optsISI.XLIM = [0,20];
optsISI.BIN_SIZE = 0.5/1000;
optsISI.DISPLAY_TEXT = 1;

ISIPlot = plotInterspikeIntervalHistogram(cds,optsExtract.NEURON_NUMBER,optsISI);




%% whole array analysis
% extract data across whole array
tic

optsExtract.STIMULI_RESPONSE = 'all';
optsExtract.STIM_ELECTRODE = unique(cds.waveforms.chanSent);
optsExtract.STIMULATIONS_PER_TRAIN = 1;
optsExtract.STIMULATION_BATCH_SIZE = 1000;

optsExtract.PRE_TIME = 20/1000; % made negative in the function
optsExtract.POST_TIME = 80/1000;
optsExtract.BIN_SIZE = 0.2/1000;
optsExtract.TIME_AFTER_STIMULATION_WAVEFORMS = 10/1000;

arrayData = extractDataAroundStimulationsWholeArray(cds,mapFileName,optsExtract);
toc


%% heatmap across whole array

% opts.STIM_ELECTRODE_PLOT = [1:numel(unique(cds.waveforms.chanSent))];
opts.STIM_ELECTRODE_PLOT = 1;
opts.WAVEFORM_TYPES_PLOT = unique(cds.waveforms.waveSent);

opts.BASELINE_PRE_TIME = -20/1000;
opts.BASELINE_POST_TIME = -3/1000;
opts.STIM_PRE_TIME = 1.2/1000;
opts.STIM_POST_TIME = 5/1000;
plotHeatmaps(arrayData,mapFileName,opts);

opts.MAX_RATIO = 1.0;
opts.MIN_RATIO = -0.05;
opts.LOG_SCALE = 1;
opts.LOG_PARAM = 4;

%% make gif
t = 1:0.1:2.1;
for i = t
    opts.STIM_PRE_TIME = i/1000;
    opts.STIM_POST_TIME = (i+1)/1000;
    opts.MAKE_BAR_PLOT = 0;
    fig = plotHeatmaps(arrayData,mapFileName,opts);
    frame = getframe(fig{1});
    im = frame2im(frame);
    [imind(:,:,1,counter),cm] = rgb2ind(im,256);
    
    
    counter = counter+1;

end

imwrite(imind,cm,'TEST_GIF.gif','gif','loopcount',inf);
%% amplitude vs. distance curve for each condition
opts.STIM_ELECTRODE_PLOT = [1:numel(unique(cds.waveforms.chanSent))];
% opts.STIM_ELECTRODE_PLOT = 1;
opts.WAVEFORM_TYPES_PLOT = unique(cds.waveforms.waveSent);

opts.BASELINE_PRE_TIME = -10/1000;
opts.BASELINE_POST_TIME = -3/1000;
opts.STIM_PRE_TIME = 1.2/1000;
opts.STIM_POST_TIME = 5/1000;

opts.LOG_SCALE = 0;
opts.PLOT_ON_ONE_FIGURE = 0;
[~,FITS]=plotAmplitudeVsDistance(arrayData,mapFileName,opts);

%% PSTH across whole array







