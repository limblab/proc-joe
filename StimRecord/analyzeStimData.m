%% set file names 

    inputData.folderpath = 'E:\Data\Joseph\Han_stim_data\Han_20191111_longTrains_dukeGen2\chan94\';
    inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

    % inputData.mapFileName = 'mapFileR:\limblab-archive\Retired Animal Logs\Monkeys\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
%     inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';


    folderpath = inputData.folderpath; % rest of code uses folderpath currently...may have switched this, not 100% certain

    inputData.task='taskCObump';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayLeftS1'; 
    inputData.monkey='monkeyHan';
    inputData.labnum = 6;

    pwd=cd;
    cd(inputData.folderpath)
    fileList = dirSorted('*spikesExtracted.nev*');
    stimInfoFileList = dirSorted('*stimInfo.mat');


%% extract relevant data for all units -- recommend saving arrayData after this step
    tic

    optsExtract.STIMULI_RESPONSE = 'all';
    optsExtract.STIMULATIONS_PER_TRAIN = 1;
    optsExtract.STIMULATION_BATCH_SIZE = 1000;
    optsExtract.DOWNSAMPLE_STIM_TIMES = 1;
    
    optsExtract.USE_STIM_CODE = 0;
    optsExtract.STIM_ELECTRODE = {};
    optsExtract.CHAN_LIST = {};

    optsExtract.PRE_TIME = 5000/1000; % made negative in the function
    optsExtract.POST_TIME = 20000/1000;

    optsExtract.BIN_SIZE = 100/1000;
    optsExtract.TIME_AFTER_STIMULATION_WAVEFORMS = 10/1000;
    optsExtract.USE_ANALOG_FOR_STIM_TIMES = 1; % this uses the analog sync line to get stim times, not sure why you would want to do anything else
    
    optsExtract.GET_KIN = 1;
    optsExtract.GET_FORCE = 0;
    
    arrayData = extractDataAroundStimulations(inputData,fileList,stimInfoFileList,optsExtract);

    toc
    
    
%% pick a unit (index in array data)
% plot raster, and PSTH for the given unit above
for arrIdx = 1:numel(arrayData)
% arrIdx = 13;
    % plot raster, and PSTH for the given unit above

%     optsPlotFunc.BIN_SIZE = optsExtract.BIN_SIZE;
    optsPlotFunc.BIN_SIZE = mode(diff(arrayData{1}.binEdges{1,1}));
    optsPlotFunc.FIGURE_SAVE = 0;
    optsPlotFunc.FIGURE_DIR = inputData.folderpath;
    optsPlotFunc.FIGURE_PREFIX = 'Han_20190923';

    optsPlotFunc.PRE_TIME = 5000/1000;
    optsPlotFunc.POST_TIME = 20000/1000;

    optsPlotFunc.SORT_DATA = '';

    optsPlotFunc.PLOT_AFTER_STIMULATION_END = 0;
    optsPlotFunc.STIMULATION_LENGTH = [];
    
%     rasterPlots = plotRasterStim(arrayData{arrIdx},arrayData{arrIdx}.NN,optsPlotFunc);

    optsPlotFunc.PLOT_ALL_ONE_FIGURE = 0;
    optsPlotFunc.PLOT_LINE = 1;
    optsPlotFunc.PLOT_TITLE = 1;    
    optsPlotFunc.PLOT_ALL_WAVES_ONE_FIGURE = 0;
% %     
    PSTHPlots = plotPSTHStim(arrayData{arrIdx},1,optsPlotFunc);

end

%% plot grid
    optsGrid.STIM_ELECTRODE = unique(arrayData{arrIdx}.CHAN_SENT);
    optsGrid.RECORDING_ELECTRODE = arrayData{arrIdx}.CHAN_REC;

    optsGrid.STIM_ELECTRODE_COLOR = {'k','r','b',[0,0.5,0],'m'};
    optsGrid.STIM_ELECTRODE_LABEL = 'string';

    optsGrid.RECORDING_ELECTRODE_COLOR = 'k';
    optsGrid.RECORDING_ELECTRODE_LABEL = 'string';

    optsGrid.FIGURE_SAVE = 0;

    optsGrid.FIGURE_DIR = '';
    optsGrid.FIGURE_PREFIX = 'Duncan_20190213_';

    ArrayPlots = plotArrayMap(arrayData{arrIdx},inputData.mapFileName(8:end),optsGrid);



%% below are whole array analyses
%% label as excitatory or inhibitory based on Hao 2016 (which is based on 
% Kraskov 2011)
% this is buggy and not really great.....probably should implement a different
% method that finds the best period of time? or deals with multiple peaks?
    optsLabel = [];
    [arrayData,report] = labelResponse(arrayData,optsLabel);


%% heatmap across whole array
disp('start')
    inputData.mapFileName = 'mapFileZ:\Basic_Sciences\Phys\L_MillerLab\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    inputData.folderpath = 'C:\Users\joh8881\Desktop\Han_20190930_trains_noAmp\';
    
    opts.STIM_ELECTRODE_PLOT = [1];
    %opts.STIM_ELECTRODE_PLOT = 1;
    %opts.WAVEFORM_TYPES_PLOT = unique(arrayData{1}.WAVEFORM_SENT);
    opts.WAVEFORM_TYPES_PLOT = [1:size(arrayData{1}.binEdges,2)];

    opts.ALL_NEURONS = 1; % 1 = plot all neurons for each stim chan, 0 = plot all stim chans for a neuron

    opts.BASELINE_PRE_TIME = -100/1000;
    opts.BASELINE_POST_TIME = -5/1000;
    opts.STIM_PRE_TIME = 0/1000;
    opts.STIM_POST_TIME = 120/1000;
    
    %time window for standardized values
    opts.PRE_STIM_WINDOW = 120/1000;
    opts.POST_STIM_WINDOW = 120/1000;

    opts.AUTO_WINDOW = 0; % 
    opts.INHIBITORY = 0;
    opts.EXCITATORY = 0;

    opts.MAX_RATIO = 1;
    opts.MIN_RATIO = -1;
    opts.LOG_SCALE = 0;
    opts.LOG_PARAM = 9;

    opts.RELATIVE_INHIBITION = 0;

    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = inputData.folderpath;
    opts.FIGURE_PREFIX = 'Han_20190924';
    
        [heatmaps, heatmap_data] = plotHeatmaps(arrayData,inputData.mapFileName(8:end),opts);
        

%% amplitude vs. distance curve for each condition -- excitation
    opts.STIM_ELECTRODE_PLOT = [1:numel(unique(arrayData{1}.CHAN_SENT))];
    opts.WAVEFORM_TYPES_PLOT = 1;%unique(arrayData{1}.WAVEFORM_SENT);
    % opts.WAVEFORM_TYPES_PLOT =;

    % maybe automatically find these
    opts.BASELINE_PRE_TIME = -15/1000;
    opts.BASELINE_POST_TIME = -5/1000;
    opts.STIM_PRE_TIME = 2/1000;
    opts.STIM_POST_TIME = 6/1000;

    opts.AUTO_WINDOW = 0; % doesn't really work...keep as 0
    opts.PLOT_SIGNIFICANT_ONLY = 0;

    opts.FIGURE_SAVE = 0;
    opts.PLOT_ON_ONE_FIGURE = 1;
    [f,FITS,data_all]=plotAmplitudeVsDistance(arrayData,inputData.mapFileName(8:end),opts);

%% latency excitation vs distance
    opts = [];
    opts.PLOT_FIT = 0;
    [peak_latency_data,~,FITS,gof] = plotLatencyExciteVsDistance(arrayData,inputData.mapFileName(8:end),opts);

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
