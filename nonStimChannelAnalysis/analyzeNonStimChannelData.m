%% set folder names 
    inputData.folderpath = 'C:\Users\joh8881\Desktop\Han_20190923_trains_noAmp\';
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
    stimInfoFileList = dirSorted('*stimInfo*');


%% extract relevant data for all units -- recommend saving arrayData after this step
    tic

    optsExtract.STIMULI_RESPONSE = 'all';
    optsExtract.STIMULATIONS_PER_TRAIN = 11;
    optsExtract.STIMULATION_BATCH_SIZE = 1000;

    optsExtract.USE_STIM_CODE = 0;
    optsExtract.STIM_ELECTRODE = {};
    optsExtract.CHAN_LIST = [];

    optsExtract.PRE_TIME = 150/1000; % made negative in the function
    optsExtract.POST_TIME = 500/1000;

    optsExtract.BIN_SIZE = 5/1000;
    optsExtract.TIME_AFTER_STIMULATION_WAVEFORMS = 10/1000;
    optsExtract.USE_ANALOG_FOR_STIM_TIMES = 1; % this uses the analog sync line to get stim times, not sure why you would want to do anything else
    optsExtract.GET_KIN = 1;
    optsExtract.GET_FORCE = 0;
    
    arrayData = extractDataAroundStimulations(inputData,fileList,stimInfoFileList,optsExtract);

    toc
    
    
%% pick a unit (index in array data)
% plot raster, and PSTH for the given unit above
for arrIdx = 6%:numel(arrayData)
% arrIdx = 1;
    % plot raster, and PSTH for the given unit above

%     optsPlotFunc.BIN_SIZE = optsExtract.BIN_SIZE;
    optsPlotFunc.BIN_SIZE = mode(diff(arrayData{1}.binEdges{1,1}));
    optsPlotFunc.FIGURE_SAVE = 0;
    optsPlotFunc.FIGURE_DIR = inputData.folderpath;
    optsPlotFunc.FIGURE_PREFIX = 'Han_20190923';

    optsPlotFunc.PRE_TIME = 50/1000;
    optsPlotFunc.POST_TIME = 200/1000;
    optsPlotFunc.SORT_DATA = '';

    optsPlotFunc.PLOT_AFTER_STIMULATION_END = 0;
    optsPlotFunc.STIMULATION_LENGTH = 0.1+0.5+0.53;
    
    rasterPlots = plotRasterStim(arrayData{arrIdx},arrayData{arrIdx}.NN,optsPlotFunc);

    optsPlotFunc.PLOT_ALL_ONE_FIGURE = 0;
    optsPlotFunc.PLOT_LINE = 1;
    optsPlotFunc.PLOT_TITLE = 1;    
    optsPlotFunc.PLOT_ALL_WAVES_ONE_FIGURE = 0;
% %     
%     PSTHPlots = plotPSTHStim(arrayData{arrIdx},1,optsPlotFunc);

end


%% heatmap across whole array

%     inputData.folderpath = 'C:\Users\joh8881\Desktop\Han_20190930_trains_noAmp\';
    inputData.folderpath = 'E:\Data\Joseph\Han_stim_data\Han_20191022_trains_noDukeAmp';
    
%    inputData.mapFileName = 'mapFileZ:\Basic_Sciences\Phys\L_MillerLab\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%    inputData.mapFileName = 'mapFileZ:\Basic_Sciences\Phys\L_MillerLab\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
   inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%    inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';

    
    opts.STIM_ELECTRODE_PLOT = [1:size(arrayData{1}.binEdges,1)];
    opts.WAVEFORM_TYPES_PLOT = [1:size(arrayData{1}.binEdges,2)];

    opts.ALL_NEURONS = 0; % 1 = plot all neurons for each stim chan, 0 = plot all stim chans for a neuron

    %time window for standardized values
    opts.BASELINE_PRE_TIME = -120/1000;
    opts.BASELINE_POST_TIME = -5/1000;
    opts.STIM_PRE_TIME = 0/1000;
    opts.STIM_POST_TIME = 120/1000;

    opts.AUTO_WINDOW = 0; % 
    opts.INHIBITORY = 0;
    opts.EXCITATORY = 0;

    opts.MAX_RATIO = 5;
    opts.MIN_RATIO = -1;
    opts.LOG_SCALE = 0;
    opts.LOG_PARAM = 9;

    opts.RELATIVE_INHIBITION = 0;

    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = inputData.folderpath;
    opts.FIGURE_PREFIX = 'Han_20191021';
        [heatmaps, heatmap_data] = plotHeatmaps(arrayData,inputData.mapFileName(8:end),opts);
        
        
        
%% compare stim response and neural response during a task
% must make sure PD related data is here, usually from running
% analyzeCObump or loading a corresponding file

%create array of averageCombinedData to pass through
%plotCombinedAverageData()

    inputData.td_all = td_all;
    inputData.mapData = mapData;
    inputData.arrayData = arrayData;
    inputData.figPrefix = 'Han_20191022';
    
    inputData.plotHeatmap = 1; %boolean 1 (plot heatmaps) or 0 (don't)
    
    test_angles = 0;
    
    averageCombinedData = nan(numel(heatmap_data),numel(test_angles));

    for i=1:numel(heatmap_data) %if multiple waveforms
        inputData.angleNumber = test_angles(ang); %if multiple 10x10 arrays
        inputData.mainChan = heatmap_data{i}.main_chan;

        inputData.dataRatioScaled = heatmap_data{i}.dataRatioScaled;
        inputData.PDscaled = PDscaled;
        inputData.waveform = heatmap_data{i}.wave; %if multiple waveforms

        combinedData = combineHeatmaps(inputData);
        %fill in averageCombinedData
        averageCombinedData(i,ang) = combinedData.average;
    end

    %plot average values vs. angle (1-360 degrees)
    % averageData = plotCombinedAverageData();

    
    
    




 