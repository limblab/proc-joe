function [ figHandle ] = plotPSTHLIB(xData,yData,optsPlotInput,optsSaveInput)
% takes in a set of data, plots that data with a bunch of parameters for
% making plots pretty and whatnot, returns figure handle
% data is a column vector

% data comes in as histogram data, yData and xData must have same shape.
% rows are data, columns are different plots.


    figHandle = '';
    %% deal with both sets of options
    optsPlot = configureOptionsPlot(optsPlotInput,xData,yData);
    optsSave = configureOptionsSave(optsSaveInput);
    
    %% check for size
    if(size(xData,2) ~= optsPlot.NUM_PLOTS || size(yData,2) ~= optsPlot.NUM_PLOTS || max(size(xData)) ~= max(size(yData)))
        warning('Aborting due to incorrectly sized inputs');
        return;
    end


    %% deal with making plot
    if(optsPlot.MAKE_FIGURE)
        figHandle = figure();
    else
        figHandle = gcf;
    end
    
    % smooth data if requested
    if(optsPlot.SMOOTH)
        [numSamples,numPlots] = size(yData);
        yDataSmoothed = zeros(numSamples,numPlots);
        
        kernelHalfLength = ceil(3*optsPlot.SMOOTH_STD_DEV/(optsPlot.BIN_SIZE));
        kernel = normpdf( -kernelHalfLength*(optsPlot.BIN_SIZE) : ...
                        optsPlot.BIN_SIZE : kernelHalfLength*(optsPlot.BIN_SIZE), ...
                        0, optsPlot.SMOOTH_STD_DEV )';
        normFactor = conv(kernel,ones(numSamples,1));
        
        % do the smoothing
        for idx = 1:numPlots
            auxSmoothed = conv(kernel,yData(:,idx))./normFactor;
            yDataSmoothed(:,idx) = auxSmoothed(kernelHalfLength+1:end-kernelHalfLength);
        end
        
        yData = yDataSmoothed;
    end
    
    
    % plot data
    for p = 1:optsPlot.NUM_PLOTS
        if(strcmpi(optsPlot.BAR_STYLE,'line')==1) % plot as line
            plot(xData(:,p),yData(:,p),'color',optsPlot.FACE_COLOR{p},'linewidth',optsPlot.LINE_WIDTH,'linestyle',optsPlot.LINE_STYLE);
        else % plot as bar
            bar(xData(:,p),yData(:,p),optsPlot.WIDTH,'EdgeColor',optsPlot.EDGE_COLOR{p},'FaceColor',optsPlot.FACE_COLOR{p},'linewidth',optsPlot.LINE_WIDTH);
        end
        hold on
    end
    
    % no recording box
    if(optsPlot.PLOT_NO_RECORDING_BOX)
        rectPos = [optsPlot.NO_RECORDING_WINDOW(1)*1000,0,diff(optsPlot.NO_RECORDING_WINDOW)*1000,optsPlot.Y_LIMITS(2)*10];
        rectangle('Position',rectPos,'FaceColor',optsPlot.NO_RECORDING_BOX_COLOR,'EdgeColor',optsPlot.NO_RECORDING_BOX_COLOR);
    end
    % deal with plot things
    if(strcmpi(optsPlot.TITLE,'')~=1)
        title(optsPlot.TITLE);
    end
    if(strcmpi(optsPlot.X_LABEL,'')~=1)
        xlabel(optsPlot.X_LABEL);
    end
    if(strcmpi(optsPlot.Y_LABEL,'')~=1)
        ylabel(optsPlot.Y_LABEL);
    end

    if(strcmpi(optsPlot.X_TICK,'')~=1)
        set(gca,'XTick',optsPlot.X_TICK);
    end
    if(strcmpi(optsPlot.Y_TICK,'')~=1)
        set(gca,'YTick',optsPlot.Y_TICK);
    end

    if(strcmpi(optsPlot.X_MINOR_TICK,'')~=1)
        set(gca,'XMinorTick',optsPlot.X_MINOR_TICK);
    end
    if(strcmpi(optsPlot.Y_MINOR_TICK,'')~=1)
        set(gca,'YMinorTick',optsPlot.Y_MINOR_TICK);
    end

    if(strcmpi(optsPlot.X_TICK_LABEL,'')~=1)
        set(gca,'XTickLabel',optsPlot.X_TICK_LABEL);
    end
    if(strcmpi(optsPlot.Y_TICK_LABEL,'')~=1)
        set(gca,'YTickLabel',optsPlot.Y_TICK_LABEL);
    end
    
    
    % legend related stuff
    if(strcmpi(optsPlot.LEGEND_STRING,'')~=1)
        l=legend(optsPlot.LEGEND_STRING);
        set(l,'box',optsPlot.LEGEND_BOX);
        set(l,'fontsize',16);
    end
    
    % do x and y limits last
    if(strcmpi(optsPlot.X_LIMITS,'')~=1)
        figHandle.CurrentAxes.XLim = optsPlot.X_LIMITS;
    end
    if(strcmpi(optsPlot.Y_LIMITS,'')~=1)
        figHandle.CurrentAxes.YLim = optsPlot.Y_LIMITS;
    end
    
    % format for lee
    formatForLee(figHandle);

    %% deal with saving plot
    if(optsSave.FIGURE_SAVE && strcmpi(optsSave.FIGURE_NAME,'')~=1 && strcmpi(optsSave.FIGURE_DIR,'')~=1)
        saveFiguresLIB(figHandle,optsSave.FIGURE_DIR,strcat(optsSave.FIGURE_NAME));
    end


end


function [optsPlot] = configureOptionsPlot(optsPlotInput,xData,yData)
    
    %% initialize possible fieldname variables
    optsPlot.MAKE_FIGURE = 1;
    optsPlot.X_LABEL = '';
    optsPlot.Y_LABEL = '';
    optsPlot.X_LIMITS = [min(xData(:,1)),max(xData(:,1))];
    optsPlot.Y_LIMITS = [min(yData(:,1)),max(yData(:,1))];
    optsPlot.X_TICK = '';
    optsPlot.Y_TICK = '';
    optsPlot.X_MINOR_TICK = '';
    optsPlot.Y_MINOR_TICK = '';
    optsPlot.X_TICK_LABEL = '';
    optsPlot.Y_TICK_LABEL = '';
    optsPlot.TITLE = '';
    optsPlot.BAR_STYLE = 'bar';
    optsPlot.EDGE_COLOR= {'k','r',[0 0.6 0],'b','m'};
    optsPlot.FACE_COLOR = {'k','r',[0 0.6 0],'b','m'};
    optsPlot.LINE_WIDTH = 1;
    optsPlot.LINE_STYLE = '-';
    optsPlot.LEGEND_STRING = '';
    optsPlot.LEGEND_BOX = 'off';
    optsPlot.NUM_PLOTS = 1;
    optsPlot.WIDTH = 0.9;
    optsPlot.BAR_COLOR = 'k';
    
    optsPlot.PLOT_NO_RECORDING_BOX = 0;
    optsPlot.NO_RECORDING_WINDOW = [0,10];
    optsPlot.NO_RECORDING_BOX_COLOR = [0,0,0];
    
    optsPlot.SMOOTH = 0;
    optsPlot.SMOOTH_STD_DEV = 1;
    
    optsPlot.BIN_SIZE = 0.2;
    %% check if in optsPlot and optsPlotInput, overwrite if so
    try
        inputFieldnames = fieldnames(optsPlotInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(optsPlot,inputFieldnames{fn}))
               optsPlot.(inputFieldnames{fn}) = optsPlotInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
    
    if(~iscell(optsPlot.EDGE_COLOR))
        optsPlot.EDGE_COLOR = {optsPlot.EDGE_COLOR};
    end
    if(~iscell(optsPlot.FACE_COLOR))
        optsPlot.FACE_COLOR = {optsPlot.FACE_COLOR};
    end
    
end

function [optsSave] = configureOptionsSave(optsSaveInput)

    optsSave.FIGURE_SAVE = 0;
    optsSave.FIGURE_DIR = '';
    optsSave.FIGURE_NAME = '';

    %% check if in optsSave and optsSaveInput, overwrite if so
    try
        inputFieldnames = fieldnames(optsSaveInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(optsSave,inputFieldnames{fn}))
               optsSave.(inputFieldnames{fn}) = optsSaveInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
end
