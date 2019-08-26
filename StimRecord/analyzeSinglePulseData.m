%% look at single pulse data. Describe activation thresholds, latency of
% different volleys, inhibition strength as a function of amplitude and
% other parameters.
% Need to create a cell array called arrayData where each entry is a recorded
% neuron. analyzeStimData does this


%% load in arrayData

%% rebin data if desired
    inputData.bin_size = 5; % in ms
    
    arrayData = rebinArrayData(arrayData,inputData);

%% plot rasters and psth for each condition and neuron

%     for arrIdx = 1:numel(arrayData)
    arrIdx = 1;
        % plot raster, and PSTH for the given unit above

    %     optsPlotFunc.BIN_SIZE = optsExtract.BIN_SIZE;
        optsPlotFunc.BIN_SIZE = mode(diff(arrayData{1}.binEdges{1,1}));
        optsPlotFunc.FIGURE_SAVE = 0;
        optsPlotFunc.FIGURE_DIR = inputData.folderpath;
        optsPlotFunc.FIGURE_PREFIX = 'Han_20190821_short';

        optsPlotFunc.PRE_TIME = 50/1000;
        optsPlotFunc.POST_TIME = 200/1000;
        optsPlotFunc.SORT_DATA = 'postStimuliTime';

        optsPlotFunc.PLOT_AFTER_STIMULATION_END = 0;
        optsPlotFunc.STIMULATION_LENGTH = 0.1+0.5+0.53;

        rasterPlots = plotRasterStim(arrayData{arrIdx},arrayData{arrIdx}.NN,optsPlotFunc);

        optsPlotFunc.PLOT_ALL_ONE_FIGURE = 0;
        optsPlotFunc.PLOT_LINE = 1;
        optsPlotFunc.PLOT_TITLE = 1;    
        optsPlotFunc.PLOT_ALL_WAVES_ONE_FIGURE = 0;
    %   
%         PSTHPlots = plotPSTHStim(arrayData{arrIdx},arrayData{arrIdx}.NN,optsPlotFunc);

%     end


%% get number of spikes (probability) as a function of amplitude and plot
% that for each neuron (on one plot probably)

    optsSpikesPlot.PRE_WINDOW = [-40,-5]/1000; % in s
    optsSpikesPlot.POST_WINDOW = [1,8]/1000; % in s
    optsSpikesPlot.PW1 = 200;
    optsSpikesPlot.PW2 = 200;
    optsSpikesPlot.POL = 0; % 0 is cathodic first
    spikesStruct = {};
    for unit = 1:numel(arrayData)
        if(unit == 1)
            optsSpikesPlot.MAKE_FIGURE = 1;
        else
            optsSpikesPlot.MAKE_FIGURE = 0;
        end
        [spikesStruct{unit},figure_handles] = plotSpikesPerAmplitude(arrayData{unit},optsSpikesPlot);
        hold on
        
        formatForLee(gcf)
        xlabel('Amplitude (\muA)');
        ylabel('Spikes per stimulation');
        set(gca,'fontsize',14)
       
    end
    
    
    
%% metric for inhibition
    optsInhibPlot = [];
    optsInhibPlot.PRE_WINDOW = [-40,-10];
    optsInhibPlot.POST_WINDOW = [0,150];
    optsInhibPlot.BIN_SIZE = 2;
    optsInhibPlot.KERNEL_LENGTH = 15;
    optsInhibPlot.BLANK_TIME = 10; % ms
    
    inhibStruct = {};
    for unit = 1:numel(arrayData)
        if(unit == 1)
            optsSpikesPlot.MAKE_FIGURE = 1;
        else
            optsSpikesPlot.MAKE_FIGURE = 0;
        end
        [inhibStruct{unit},figure_handles] = plotInhibitionDuration(arrayData{unit},optsInhibPlot);
%         hold on
%         
%         formatForLee(gcf)
%         xlabel('Amplitude (\muA)');
%         ylabel('Inhibition duration');
%         set(gca,'fontsize',14)
       
    end
    
    
    
    
    
    
    
    
    
    
    
    
    















