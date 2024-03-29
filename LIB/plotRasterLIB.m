function [ figHandle ] = plotRasterLIB(xData,yData,optsPlotInput,optsSaveInput)
% takes in a set of data, plots that data with a bunch of parameters for
% making plots pretty and whatnot, returns figure handle
% xdata is the ts, ydata is the row number, both should be column vectors

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
    %% sort data if requested
    if(strcmpi(optsPlot.SORT_DATA,'no')~=1)
        switch optsPlot.SORT_DATA
            case 'postStimuliTime'
                yDataShort = yData(xData > 0);
                xDataShort = xData(xData > 0);
%                 xDataShort(xData <= 0) = 100000; % make the times below 0 a big number to sort post stimuli time
%                 yDataShort = yData;
                [~,sortIdx] = sort(xDataShort);
                yDataShort = yDataShort(sortIdx);
                sortIdx = unique(yDataShort,'stable');
                % add trials with no data to end of sortIdx
                allIdx = 1:max(yData);
                sortIdx = [sortIdx; allIdx(setdiff(allIdx,sortIdx))'];
                yDataSorted = yData;
                xDataSorted = xData;
                counter = 1;
                for sIdx = 1:numel(sortIdx)
                    yDataSorted(counter:counter+sum(yData==sortIdx(sIdx))-1) = sIdx;
                    xDataSorted(counter:counter+sum(yData==sortIdx(sIdx))-1) = xData(yData==sortIdx(sIdx));
                    counter = counter+sum(yData==sortIdx(sIdx));
                end

                yData = yDataSorted;
                xData = xDataSorted;
        end

    end
    
    if(optsPlot.MAKE_FIGURE)
        figHandle = figure('Position',[680 558 620 420]);
    else
        figHandle = gcf;
    end
    if(optsPlot.MAKE_SUBPLOT)
        subplot(optsPlot.SUBPLOT_SIZE(1),optsPlot.SUBPLOT_SIZE(2),optsPlot.SUBPLOT_IDX)
    end
    % plot data
    if(strcmpi(optsPlot.MARKER_STYLE,'line')==0) % normal plot routine
        plot(xData,yData,'.','marker',optsPlot.MARKER_STYLE,'color',optsPlot.MARKER_COLOR,'markersize',optsPlot.MARKER_SIZE);
    else % plot as vertical lines
        plot([xData.';xData.';nan(1,length(xData))],...
            [(yData-optsPlot.LINE_LENGTH/2).';(yData+optsPlot.LINE_LENGTH/2).';nan(1,length(yData))],...
            'color',optsPlot.MARKER_COLOR,'linewidth',optsPlot.LINE_WIDTH)
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
            plot([-10000,10000],[yVal, yVal],'-','Color',optsPlot.DIVIDING_LINES_COLORS{idx},'linewidth',1);
        end
        
        if(strcmpi(optsPlot.DUKE_BLANKED_REGION,'')~=1)
            for idx = 1:size(optsPlot.DUKE_BLANKED_REGION,1)
                if(idx == 1)
                    y=optsPlot.Y_LIMITS(1);
                    h=optsPlot.DIVIDING_LINES(idx);
                elseif(idx == size(optsPlot.DUKE_BLANKED_REGION,1))
                    y=optsPlot.DIVIDING_LINES(idx-1);
                    h=optsPlot.Y_LIMITS(2)-optsPlot.DIVIDING_LINES(idx-1);
                else
                    y=optsPlot.DIVIDING_LINES(idx-1);
                    h=optsPlot.DIVIDING_LINES(idx)-optsPlot.DIVIDING_LINES(idx-1);
                end
                pos = [optsPlot.DUKE_BLANKED_REGION(idx,1),y,...
                    optsPlot.DUKE_BLANKED_REGION(idx,2)-optsPlot.DUKE_BLANKED_REGION(idx,1),h]; % [x,y,w,h]
                temp = rectangle('Position',pos,'edgecolor','none','facecolor',[1,0,0,0.30]);
                uistack(temp,'bottom');
            end
        end
        
        if(strcmpi(optsPlot.BLACK_BLANKED_REGION,'')~=1)
            for idx = 1:size(optsPlot.BLACK_BLANKED_REGION,1)
                if(idx == 1)
                    y=optsPlot.Y_LIMITS(1);
                    h=optsPlot.DIVIDING_LINES(idx);
                elseif(idx == size(optsPlot.BLACK_BLANKED_REGION,1))
                    y=optsPlot.DIVIDING_LINES(idx-1);
                    h=optsPlot.Y_LIMITS(2)-optsPlot.DIVIDING_LINES(idx-1);
                else
                    y=optsPlot.DIVIDING_LINES(idx-1);
                    h=optsPlot.DIVIDING_LINES(idx)-optsPlot.DIVIDING_LINES(idx-1);
                end
                pos = [optsPlot.BLACK_BLANKED_REGION(idx,1),y,...
                    optsPlot.BLACK_BLANKED_REGION(idx,2)-optsPlot.BLACK_BLANKED_REGION(idx,1),h]; % [x,y,w,h]
                temp = rectangle('Position',pos,'edgecolor','none','facecolor',[0,0,0,0.2]);
                uistack(temp,'bottom');
            end
        end
        
    end
    
    % plot stim times if prompted
    if(optsPlot.PLOT_STIM_TIME == 1 && ~isempty(optsPlot.STIM_DATA_X) && ~isempty(optsPlot.STIM_DATA_Y))
        if(strcmpi(optsPlot.MARKER_STYLE,'line')==0) % normal plot routine
            plot(optsPlot.STIM_DATA_X,optsPlot.STIM_DATA_Y,'.','marker',optsPlot.MARKER_STYLE,'color',optsPlot.STIM_DATA_COLOR,'markersize',optsPlot.MARKER_SIZE);
        else % plot as vertical lines
            plot([optsPlot.STIM_DATA_X.';optsPlot.STIM_DATA_X.';nan(1,length(optsPlot.STIM_DATA_X))],...
                [(optsPlot.STIM_DATA_Y-optsPlot.LINE_LENGTH/2).';(optsPlot.STIM_DATA_Y+optsPlot.LINE_LENGTH/2).';nan(1,length(optsPlot.STIM_DATA_Y))],...
                'color',optsPlot.STIM_DATA_COLOR,'linewidth',optsPlot.STIM_LINE_WIDTH)
        end
    end
    

    
    % do x and y limits last
    if(strcmpi(optsPlot.X_LIMITS,'')~=1)
        figHandle.CurrentAxes.XLim = optsPlot.X_LIMITS;
    end
    if(strcmpi(optsPlot.Y_LIMITS,'')~=1)
        figHandle.CurrentAxes.YLim = optsPlot.Y_LIMITS;
    end
    
    % plot zero line if prompted
    if(optsPlot.PLOT_ZERO_LINE)
        hold on
        plot([0,0],figHandle.CurrentAxes.YLim,'r--','linewidth',1)
    end
    
    % format for lee
    formatForLee(figHandle);

    set(gca,'FontSize',optsPlot.FONT_SIZE);
    
    if(strcmpi(optsPlot.X_TICK,'')~=1)
        set(gca,'XTick',optsPlot.X_TICK);
    end
    if(strcmpi(optsPlot.Y_TICK,'')~=1)
        try
            set(gca,'YTick',optsPlot.Y_TICK);
        catch
        end
    end

    if(strcmpi(optsPlot.X_MINOR_TICK,'')~=1)
        set(gca,'XMinorTick',optsPlot.X_MINOR_TICK);
    end
    if(strcmpi(optsPlot.Y_MINOR_TICK,'')~=1)
        set(gca,'YMinorTick',optsPlot.Y_MINOR_TICK);
    end
    
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
    optsPlot.LINE_LENGTH = 1;
    optsPlot.TITLE = '';
    optsPlot.LINE_WIDTH = 2;
    optsPlot.MARKER_STYLE = '.';
    optsPlot.MARKER_COLOR = 'k';
    optsPlot.MARKER_SIZE = 3;
    optsPlot.DIVIDING_LINES = '';
    optsPlot.DIVIDING_LINES_COLORS = {'k','k','k','k','k','k','k','k','k','k','k','k','k','k','k'};
    optsPlot.DUKE_BLANKED_REGION = '';
    optsPlot.BLACK_BLANKED_REGION = '';
    optsPlot.PLOT_STIM_TIME = 0;
    optsPlot.STIM_DATA_X = [];
    optsPlot.STIM_DATA_Y = [];
    optsPlot.STIM_DATA_COLOR = 'r';
    optsPlot.STIM_LINE_WIDTH = 1.2;
    optsPlot.SORT_DATA = 'no'; % postStimuliTime
    optsPlot.FONT_SIZE = 14;
    
    optsPlot.MAKE_SUBPLOT = 0;
    optsPlot.SUBPLOT_IDX = -1;
    optsPlot.SUBPLOT_SIZE = [];
    optsPlot.PLOT_ZERO_LINE = 0;
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