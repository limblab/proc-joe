function [figureHandle] = plotPSTHStim(unitData,NEURON_NUMBER,optsPlot)



    %% configure opts and set default values
    opts = configureOpts(optsPlot);
    
    NUM_WAVEFORM_TYPES = size(unitData.spikeTrialTimes,2);
    NUM_CHANS = size(unitData.spikeTrialTimes,1);
    CHAN_LIST = unitData.CHAN_LIST;
    
    %% format data and call general plot function
    
    %% plot data - stimuli data is y-axis. spike data is x-axis

    if(opts.PLOT_ALL_ONE_FIGURE)
        NUM_CHANS_PLOT = 1;
        NUM_WAVEFORM_TYPES_PLOT = 1;
    elseif(opts.PLOT_ALL_WAVES_ONE_FIGURE)
        NUM_CHANS_PLOT = NUM_CHANS;
        NUM_WAVEFORM_TYPES_PLOT = 1;
    else
        NUM_CHANS_PLOT = NUM_CHANS;
        NUM_WAVEFORM_TYPES_PLOT = NUM_WAVEFORM_TYPES;
    end
    
    for chan = 1:NUM_CHANS_PLOT
        for wave = 1:NUM_WAVEFORM_TYPES_PLOT

            % get data
            xData = [];
            yData = [];
            if(opts.PLOT_ALL_ONE_FIGURE)
                optsPlot.NUM_PLOTS = NUM_CHANS*NUM_WAVEFORM_TYPES;
                for chanTemp = 1:NUM_CHANS
                    for waveTemp = 1:NUM_WAVEFORM_TYPES
                        bEtemp = unitData.bE{chanTemp,waveTemp};
                        xDataTemp = bEtemp(1:end-1)+(bEtemp(2)-bEtemp(1))/2;
                        xDataTemp = xDataTemp';
                        yDataTemp = unitData.bC{chanTemp,waveTemp};
                        yDataTemp = yDataTemp';
                        xData(:,end+1) = xDataTemp;
                        yData(:,end+1) = yDataTemp;
                    end
                end
            elseif(opts.PLOT_ALL_WAVES_ONE_FIGURE)
                optsPlot.NUM_PLOTS = NUM_WAVEFORM_TYPES;
                for waveTemp = 1:NUM_WAVEFORM_TYPES
                    bEtemp = unitData.bE{chan,waveTemp};
                    xDataTemp = bEtemp(1:end-1)+(bEtemp(2)-bEtemp(1))/2;
                    xDataTemp = xDataTemp';
                    yDataTemp = unitData.bC{chan,waveTemp};
                    yDataTemp = yDataTemp';
                    xData(:,end+1) = xDataTemp;
                    yData(:,end+1) = yDataTemp;
                end
            else
                optsPlot.NUM_PLOTS = 1;

                bEtemp = unitData.bE{chan,wave};
                xData = bEtemp(1:end-1)+(bEtemp(2)-bEtemp(1))/2;
                xData = xData';
                yData = unitData.bC{chan,wave};
                yData = yData';
            end
            
            optsPlot.SMOOTH = opts.SMOOTH;
            optsPlot.SMOOTH_STD_DEV = opts.SMOOTH_STD_DEV;
            optsPlot.BIN_SIZE = opts.BIN_SIZE*1000;
            
            optsPlot.FACE_COLOR = opts.LINE_COLOR;
            optsPlot.EDGE_COLOR = opts.LINE_COLOR;
            
            % format for plotting, then plot
            if(opts.PLOT_LINE)
                optsPlot.BAR_STYLE = 'line';
            else
                optsPlot.BAR_STYLE = 'bar';
            end
            if(opts.MAKE_LEGEND && ~isempty(opts.LEGEND_STR))
                optsPlot.LEGEND_STRING = opts.LEGEND_STR;
                optsPlot.LEGEND_BOX = 'off';
            end      
            optsPlot.Y_LIMITS = [0,unitData.binMaxYLim];
            optsPlot.X_LIMITS = [-opts.PRE_TIME*1000,opts.POST_TIME*1000];
            optsPlot.Y_LABEL = 'Number of spikes per stimulation';
            if(~opts.MAKE_SUBPLOTS || wave == NUM_WAVEFORM_TYPES)
                optsPlot.X_LABEL = 'Time after stimulation onset (ms)';
            end
            if(opts.PLOT_TITLE)
                if(strcmp(opts.TITLE_TO_PLOT,'') == 0)
                    optsPlot.TITLE = titleToPlot;
                elseif(opts.PLOT_ALL_ONE_FIGURE)
                    % no title
                else
                    optsPlot.TITLE = strcat('Stim Chan: ',num2str(CHAN_LIST{chan}),' Wave: ',num2str(wave));
                end
            end

            optsPlot.MAKE_FIGURE = opts.MAKE_FIGURE;

            optsSave.FIGURE_SAVE = opts.FIGURE_SAVE;
            optsSave.FIGURE_DIR = opts.FIGURE_DIR;
            
            if(~opts.MAKE_SUBPLOTS && opts.FIGURE_SAVE && opts.PLOT_ALL_ONE_FIGURE)
                optsSave.FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'nn',num2str(NEURON_NUMBER),'_chan',num2str(unitData.CHAN_REC),'_ALL_PSTH');
            elseif(~opts.MAKE_SUBPLOTS && opts.FIGURE_SAVE)
                optsSave.FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'nn',num2str(NEURON_NUMBER),'_chan',num2str(unitData.CHAN_REC),'_stimChan',num2str(unitData.CHAN_LIST{chan}),'_waveNum',num2str(wave),'_PSTH');
            end
            
            optsPlot.PLOT_NO_RECORDING_BOX = opts.PLOT_NO_RECORDING_BOX;
            optsPlot.NO_RECORDING_WINDOW = opts.NO_RECORDING_WINDOW;
            optsPlot.NO_RECORDING_BOX_COLOR = opts.NO_RECORDING_BOX_COLOR;
            
            figureHandle{chan,wave} = plotPSTHLIB(xData,yData,optsPlot,optsSave);


        end

    end
    




end


function [opts] = configureOpts(optsInput)

    opts.PLOT_ALL_ONE_FIGURE = 0;
    opts.PLOT_ALL_WAVES_ONE_FIGURE = 0;
    
    opts.PLOT_TITLE = 1;
    opts.TITLE_TO_PLOT = '';
    
    opts.SMOOTH = 0;
    opts.SMOOTH_STD_DEV = 1;
    
    opts.PLOT_LINE = 0;
    opts.LINE_COLOR = {getColorFromList(1,1),getColorFromList(1,2),getColorFromList(1,3),getColorFromList(1,4),getColorFromList(1,5),getColorFromList(1,6)};
    opts.PLOT_NO_RECORDING_BOX = 0;
    opts.NO_RECORDING_WINDOW = [0,1.0/1000];
    opts.NO_RECORDING_BOX_COLOR = [1.0,0.6,0.6];
    
    opts.MAKE_LEGEND = 0;
    opts.LEGEND_STR = '';
    
    opts.MAKE_FIGURE = 1;
    opts.MAKE_SUBPLOTS = 0;
    
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = '';
    opts.FIGURE_PREFIX = '';
    
    opts.PRE_TIME = 5/1000;
    opts.POST_TIME = 30/1000;
    opts.BIN_SIZE = 0.2/1000;
    
    %% check if in optsSave and optsSaveInput, overwrite if so
    try
        inputFieldnames = fieldnames(optsInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(opts,inputFieldnames{fn}))
               opts.(inputFieldnames{fn}) = optsInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
end