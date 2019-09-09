%% look at single pulse data. Describe activation thresholds, latency of
% different volleys, inhibition strength as a function of amplitude and
% other parameters.
% Need to create a cell array called arrayData where each entry is a recorded
% neuron. analyzeStimData does this


%% load in arrayData

%% rebin data if desired
    inputData.bin_size = 1; % in ms
    
    arrayData = rebinArrayData(arrayData,inputData.bin_size);

%% plot rasters and psth for each condition and neuron

%     for arrIdx = 1:numel(arrayData)
    arrIdx = 2;
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
    optsSpikesPlot.POST_WINDOW = [1,10]/1000; % in s
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
    optsInhibPlot.PRE_WINDOW = [-100,-10];
    optsInhibPlot.POST_WINDOW = [0,200];
    optsInhibPlot.MAX_TIME_START = 40; % ms
    optsInhibPlot.BIN_SIZE = 1;
    optsInhibPlot.KERNEL_LENGTH = 5;
    optsInhibPlot.BLANK_TIME = 5; % ms
    
    optsInhibPlot.PW1 = 200;
    optsInhibPlot.PW2 = 200;
    optsInhibPlot.POL = 0; % 0 is cathodic first
    
    
    inhibStruct = {};
    figure();
    for unit = 1:numel(arrayData)
        if(unit == 1)
            optsSpikesPlot.MAKE_FIGURE = 1;
        else
            optsSpikesPlot.MAKE_FIGURE = 0;
        end
        [inhibStruct{unit},figure_handles] = plotInhibitionDuration(arrayData{unit},optsInhibPlot);
        hold on
        formatForLee(gcf)
        xlabel('Amplitude (\muA)');
        ylabel('Inhibition duration (ms)');
        set(gca,'fontsize',14)
       
    end
    
    
%% linear model predicting inhib duration based on amp and unit
    amp = [];
    unit = [];
    inhib_dur = [];
    for u = 1:numel(inhibStruct)
        amp = [amp,inhibStruct{u}.amp];
        inhib_dur = [inhib_dur,inhibStruct{u}.inhib_dur'];
        unit = [unit,u*ones(size(inhibStruct{u}.inhib_dur))'];
        
    end
    
    % remove nan's
    keep_mask = ~isnan(amp) & ~isnan(inhib_dur) & ~isnan(unit);
    amp = amp(keep_mask);
    inhib_dur = inhib_dur(keep_mask);
    unit = unit(keep_mask);
    
    tbl = table(amp',inhib_dur',unit','VariableNames',{'amp','inhib_dur','unit'});
    mdl = fitlm(tbl,'inhib_dur~amp+unit')
    
    
%%
unit_idx = 2; cond = 2;
    figure();
    plot(inhibStruct{unit_idx}.filtered_PSTH(cond,:))
    disp(inhibStruct{unit_idx}.threshold(cond))
    
    
    
    
    
    
    
    















