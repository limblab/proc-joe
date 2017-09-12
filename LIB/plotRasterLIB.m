function [ figHandle ] = plotRasterLIB(xData,yData,optsPlotInput,optsSaveInput)
% takes in a set of data, plots that data with a bunch of parameters for
% making plots pretty and whatnot, returns figure handle


    %% deal with both sets of options

    optsPlot = configureOptionsPlot(optsPlotInput,xData,yData);
    optsSave = configureOptionsSave(optsSaveInput);
    
    figHandle = '';
    %% reshape data to column vectors
    if(min(size(xData)) > 1 || min(size(yData))>1 || max(size(xData)) ~= max(size(yData)))
        warning('Aborting due to incorrectly sized inputs');
        return;
    end
    if(~iscolumn(xData))
        xData = xData';
    end
    if(~iscolumn(yData))
        yData = yData';
    end

    %% deal with making plot
    if(optsPlot.MAKE_FIGURE)
        figHandle = figure();
    else
        figHandle = gcf;
    end
    % plot data
    if(strcmpi(optsPlot.MARKER_STYLE,'line')==0) % normal plot routine
        plot(xData,yData,'.','marker',optsPlot.MARKER_STYLE,'color',optsPlot.MARKER_COLOR,'markersize',optsPlot.MARKER_SIZE);
    else % plot as vertical lines
        plot([xData,xData],...
            [yData-optsPlot.LINE_LENGTH/2,yData+optsPlot.LINE_LENGTH/2]',...
            optsPlot.MARKER_COLOR,'linewidth',optsPlot.LINE_WIDTH)
    end

    % deal with plot things
    if(strcmpi(optsPlot.X_LIMITS,'')~=1)
        figHandle.CurrentAxes.XLim = optsPlot.X_LIMITS;
    end
    if(strcmpi(optsPlot.Y_LIMITS,'')~=1)
        figHandle.CurrentAxes.YLim = optsPlot.Y_LIMITS;
    end
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
    
    % dividing lines if prompted
    if(strcmpi(optsPlot.DIVIDING_LINES,'')~=1)
        for idx = 1:numel(optsPlot.DIVIDING_LINES)
            hold on
            yVal = optsPlot.DIVIDING_LINES(idx);
            plot([figHandle.CurrentAxes.XLim],[yVal, yVal],'-','Color',optsPlot.DIVIDING_LINES_COLORS{idx});
        end
    end
    % format for lee
    formatForLee(figHandle);

    %% deal with saving plot
    if(optsSave.FIGURE_SAVE && strcmpi(optsSave.FIGURE_NAME,'')~=1 && strcmpi(optsSave.FIGURE_DIR,'')~=1)
        saveFiguresLIB(figHandle,optsSave.FIGURE_DIR,optsSave.FIGURE_NAME);
    end




end


function [optsPlot] = configureOptionsPlot(optsPlotInput,xData,yData)
    
    %% initialize possible fieldname variables
    optsPlot.MAKE_FIGURE = 1;
    optsPlot.X_LABEL = '';
    optsPlot.Y_LABEL = '';
    optsPlot.X_LIMITS = [min(xData),max(xData)];
    optsPlot.Y_LIMITS = [min(yData),max(yData)];
    optsPlot.X_TICK = '';
    optsPlot.Y_TICK = '';
    optsPlot.X_MINOR_TICK = '';
    optsPlot.Y_MINOR_TICK = '';
    optsPlot.X_TICK_LABEL = '';
    optsPlot.Y_TICK_LABEL = '';
    optsPlot.LINE_STYLE = '';
    optsPlot.TITLE = '';
    optsPlot.MARKER_STYLE = '.';
    optsPlot.MARKER_COLOR = 'k';
    optsPlot.MARKER_SIZE = 4;
    optsPlot.DIVIDING_LINES = '';
    optsPlot.DIVIDING_LINES_COLORS = '';
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
     
    if(numel(optsPlot.DIVIDING_LINES) > numel(optsPlot.DIVIDING_LINES_COLORS))
        diff = numel(optsPlot.DIVIDING_LINES) - optsPlot.DIVIDING_LINES_COLORS;
        for d = 1:diff
            optsPlot.DIVIDING_LINES_COLORS{end+1,1} = 'k';
        end
    elseif(iscell(optsPlot.DIVIDING_LINES) && ~iscell(optsPlot.DIVIDING_LINES_COLORS,'')==1)
        diff = numel(optsPlot.DIVIDING_LINES);
        for d = 1:diff
            optsPlot.DIVIDING_LINES_COLORS{end+1,1} = 'k';
        end
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