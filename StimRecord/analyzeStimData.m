%% set file names 
inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Chips_20171026\processed\';
inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
folderpath = inputData.folderpath; % rest of code uses folderpath currently

inputData.task='taskCObump';
inputData.ranBy='ranByJoseph'; 
inputData.array1='arrayLeftS1'; 
inputData.monkey='monkeyChips';
inputData.labnum = 6;

pwd=cd;
cd(inputData.folderpath)
fileList = dirSorted('*spikesExtracted.nev*');
stimInfoFileList = dirSorted('*stimInfo*');


%% extract relevant data for all units -- highly recommend saving arrayData after this step
tic

optsExtract.STIMULI_RESPONSE = 'all';
optsExtract.STIMULATIONS_PER_TRAIN = 1;
optsExtract.STIMULATION_BATCH_SIZE = 1000;

optsExtract.PRE_TIME = 20/1000; % made negative in the function
optsExtract.POST_TIME = 80/1000;
optsExtract.BIN_SIZE = 0.2/1000;
optsExtract.TIME_AFTER_STIMULATION_WAVEFORMS = 10/1000;

optsExtract.USE_ANALOG_FOR_STIM_TIMES = 1;

arrayData = extractDataAroundStimulations(inputData,fileList,stimInfoFileList,optsExtract);

toc
%% pick a unit (index in array data)
arrIdx = 18;

% plot raster, and PSTH for the given unit above

optsPlotFunc.BIN_SIZE = optsExtract.BIN_SIZE;
optsPlotFunc.FIGURE_SAVE = 0;
optsPlotFunc.FIGURE_DIR = folderpath;
optsPlotFunc.FIGURE_PREFIX = 'Chips_20171024_short_';

optsPlotFunc.PRE_TIME = 15/1000;
optsPlotFunc.POST_TIME = 60/1000;
optsPlotFunc.SORT_DATA = 'postStimuliTime';

rasterPlots = plotRasterStim(arrayData{arrIdx},arrayData{arrIdx}.NN,optsPlotFunc);

optsPlotFunc.PLOT_ALL_ONE_FIGURE = 1;
optsPlotFunc.PLOT_LINE = 1;
optsPlotFunc.PLOT_TITLE = 0;

PSTHPlots = plotPSTHStim(arrayData{arrIdx},arrayData{arrIdx}.NN,optsPlotFunc);


%% plot grid
optsGrid.STIM_ELECTRODE = unique(arrayData{arrIdx}.CHAN_SENT);
optsGrid.RECORDING_ELECTRODE = arrayData{arrIdx}.CHAN_REC;

optsGrid.STIM_ELECTRODE_COLOR = {'k','r','b',[0,0.5,0],'m'};
optsGrid.STIM_ELECTRODE_LABEL = 'string';

optsGrid.RECORDING_ELECTRODE_COLOR = 'k';
optsGrid.RECORDING_ELECTRODE_LABEL = 'string';

optsGrid.FIGURE_SAVE = 0;

optsGrid.FIGURE_DIR = '';
optsGrid.FIGURE_PREFIX = 'Chips_20171026_';

ArrayPlots = plotArrayMap(arrayData{arrIdx},inputData.mapFileName(8:end),optsGrid);



%% below are whole array analyses
%% label as excitatory or inhibitory based on Hao 2016 (which is based on 
% Kraskov 2011)
%% this is buggy and not really great.....probably should implement a different
% method that finds the best period of time? or deals with multiple peaks?
optsLabel.BRUTE_FORCE_MAX_WINDOW = 0;
optsLabel.MIN_WINDOW_SIZE = 3;
optsLabel.BRUTE_FORCE_EXCITE_GAMMA = 3;
optsLabel.BRUTE_FORCE_INHIB_GAMMA = 0.5;

[arrayData,report] = labelResponse(arrayData,optsLabel);


%% heatmap across whole array

opts.STIM_ELECTRODE_PLOT = [1:numel(unique(arrayData{1}.CHAN_SENT))];
% opts.STIM_ELECTRODE_PLOT = 1;
opts.WAVEFORM_TYPES_PLOT = unique(arrayData{1}.WAVEFORM_SENT);

opts.BASELINE_PRE_TIME = -20/1000;
opts.BASELINE_POST_TIME = -2/1000;
opts.STIM_PRE_TIME = 1/1000;
opts.STIM_POST_TIME = 10/1000;

opts.AUTO_WINDOW = 0; % 
opts.INHIBITORY = 1;
opts.EXCITATORY = 0;

opts.MAX_RATIO = 1.0;
opts.MIN_RATIO = -0.1;
opts.LOG_SCALE = 1;
opts.LOG_PARAM = 9;

opts.RELATIVE_INHIBITION = 0;

opts.FIGURE_SAVE = 0;
opts.FIGURE_DIR = '';
opts.FIGURE_PREFIX = 'Chips_20171026';
    plotHeatmaps(arrayData,inputData.mapFileName(8:end),opts);

%% amplitude vs. distance curve for each condition -- excitation
opts.STIM_ELECTRODE_PLOT = [1:numel(unique(arrayData{1}.CHAN_SENT))];
opts.WAVEFORM_TYPES_PLOT = unique(arrayData{1}.WAVEFORM_SENT);

% maybe automatically find these
opts.BASELINE_PRE_TIME = -20/1000;
opts.BASELINE_POST_TIME = -3/1000;
opts.STIM_PRE_TIME = 1.2/1000;
opts.STIM_POST_TIME = 5/1000;

opts.AUTO_WINDOW = 1;

opts.FIGURE_SAVE = 0;
opts.PLOT_ON_ONE_FIGURE = 1;
[~,FITS]=plotAmplitudeVsDistance(arrayData,inputData.mapFileName(8:end),opts);

%% latency excitation vs distance
[~,FITS] = plotLatencyExciteVsDistance(arrayData,inputData.mapFileName(8:end),opts);

%% inhibition duration vs distance
[~,FITS] = plotInhibitionDurationVsDistance(arrayData,inputData.mapFileName(8:end),opts);

%% plot amplitude of response vs. inhibition duration
% opts.UNITS = [5,6,7,10,11,13,14,15,19,20];
% opts.UNITS = [];
[~,FITS] = plotAmplitudeVsInhibitionDuration(arrayData,inputData.mapFileName(8:end),opts);


%% make gif
t = 1:0.5:50;
counter = 1;
imind = [];
opts.STIM_ELECTRODE_PLOT = 1;
opts.WAVEFORM_TYPES_PLOT = 3;

opts.BASELINE_PRE_TIME = -10/1000;
opts.BASELINE_POST_TIME = -3/1000;

for i = t
    opts.STIM_PRE_TIME = i/1000;
    opts.STIM_POST_TIME = (i+5)/1000;
    opts.MAKE_BAR_PLOT = 0;
    fig = plotHeatmaps(arrayData,mapFileName,opts);
    frame = getframe(fig{1});
    im = frame2im(frame);
    imind(:,:,:,counter) = double(im)/256.0;
    
    
    counter = counter+1;
    close all
end

v = VideoWriter('C:\Users\Joseph\Desktop\Chips_20171026_noRelativeInhibition.mp4','MPEG-4');
v.FrameRate = 30;
open(v)
writeVideo(v,imind)
close(v)



%% BROKEN -- don't currently store true spike times for every stimulation, makes it hard to compute an ISI
% fix -- either store true spike times or take each stimulation
% independently
%%  ISI
optsISI.XLIM = [0,20];
optsISI.BIN_SIZE = 0.2/1000;
optsISI.DISPLAY_TEXT = 1;

optsISI.FIGURE_SAVE = 0;
optsISI.FIGURE_PREFIX = 'Chips_20171026';
optsISI.FIGURE_DIR = folderpath;

[ISIPlot,ISIdata] = plotInterspikeIntervalHistogram(cds,optsExtract.NEURON_NUMBER,optsISI);




%% Code below requires outputData

%% plot waveforms (raw and filtered)
optsWaveforms = [];

optsWaveforms.FIGURE_SAVE = 0;

optsWaveforms.FIGURE_DIR = '';
optsWaveforms.FIGURE_PREFIX = 'Chips_20171024_';

optsWaveforms.YLIM = [-800,800];
optsWaveforms.TIME_AFTER_STIMULATION_ARTIFACT = 4/1000;
%% BROKEN WITH NEW CODE -- ALL BELOW ARE PROBABLY ALSO BROKEN
WaveformPlots = plotWaveformsStim(cds,stimInfo,rawData,optsExtract.NEURON_NUMBER,optsWaveforms);

%% plot artifacts (raw and filtered)

optsArtifacts = [];
optsArtifacts.WAVEFORMS_TO_PLOT = [5];
optsArtifacts.CHANS_TO_PLOT = [1];
optsArtifacts.MAX_WAVES_PLOT = 5;
optsArtifacts.TIME_AFTER_STIMULATION_ARTIFACT = 6/1000;
optsArtifacts.YLIM = [-600,600];
optsArtifacts.XLIM = [0,10];
optsArtifacts.ADJUST_FOR_BIT_ERROR = 0;
optsArtifacts.ROW_SUBPLOT = 2;
optsArtifacts.COL_SUBPLOT = 2;
optsArtifacts.NUM_PAD = 200;
optsArtifacts.RANDOM_SAMPLE = 1;
optsArtifacts.ARTIFACT_MULTIPLIER = 5;

optsArtifacts.PLOT_FILTERED = 1;

ArtifactPlots = plotArtifactsStim(cds,stimInfo,artifactData,optsExtract.NEURON_NUMBER,optsArtifacts);



% legacy code below
% % % % %% extract unit data for whole array -- recommend saving after this to avoid long loading times
% % % % % extract data across whole array
% % % % tic
% % % % 
% % % % optsExtract.STIMULI_RESPONSE = 'all';
% % % % optsExtract.STIM_ELECTRODE = unique(stimInfo.waveforms.chanSent);
% % % % optsExtract.STIMULATIONS_PER_TRAIN = 1;
% % % % optsExtract.STIMULATION_BATCH_SIZE = 1000;
% % % % 
% % % % optsExtract.PRE_TIME = 20/1000; % made negative in the function
% % % % optsExtract.POST_TIME = 80/1000;
% % % % optsExtract.BIN_SIZE = 0.5/1000;
% % % % optsExtract.TIME_AFTER_STIMULATION_WAVEFORMS = 10/1000;
% % % % 
% % % % arrayData = extractDataAroundStimulationsWholeArray(cds,stimInfo,mapFileName,optsExtract);
% % % % toc