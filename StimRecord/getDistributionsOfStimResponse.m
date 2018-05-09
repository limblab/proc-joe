function [arrayDataFits] = getDistributionsOfStimResponse( arrayData,opts )
% 

    opts = configureOpts(opts);
    globalXLimits = [intmax,intmin];
    globalYLimits = globalXLimits;
    
    %% get data for each window
    means = zeros(size(opts.WINDOW,1),size(arrayData,2),size(arrayData{1}.bC,1),size(arrayData{1}.bC,2));
    vars = zeros(size(opts.WINDOW,1),size(arrayData,2),size(arrayData{1}.bC,1),size(arrayData{1}.bC,2));
    
    for windowIdx = 1:size(opts.WINDOW,1) % for each window
        for arrayIdx = 1:size(arrayData,2) % for each unit
            for chan = 1:size(arrayData{arrayIdx}.bC,1) % for each channel stimulated
                for waveNum = 1:size(arrayData{arrayIdx}.bC,2) % for each waveform stimulated with
                    binIdx = [find(abs(arrayData{arrayIdx}.bE{chan,waveNum}-opts.WINDOW(windowIdx,1))<1E-4), find(abs(arrayData{arrayIdx}.bE{chan,waveNum}-opts.WINDOW(windowIdx,2))<1E-4)];
                    means(windowIdx,arrayIdx,chan,waveNum) = mean(arrayData{arrayIdx}.bC{chan,waveNum}(binIdx(1):binIdx(2)));
                    vars(windowIdx,arrayIdx,chan,waveNum) = var(arrayData{arrayIdx}.bC{chan,waveNum}(binIdx(1):binIdx(2)));
                end % end for wave
            end % end for chan
        end % end for unit
        % reorganize data to compress across all units, channels, and
        % waveforms sent
        meanWindow{windowIdx} = reshape(means(windowIdx,:,:,:),numel(means(windowIdx,:,:,:)),1);
        varWindow{windowIdx} = reshape(vars(windowIdx,:,:,:),numel(vars(windowIdx,:,:,:)),1);
        
    end % end for window
    
    for windowIdx = 1:numel(meanWindow)

        %% fit data with line going through (0,0)
        [F{windowIdx},gof{windowIdx}] = fit(meanWindow{windowIdx},varWindow{windowIdx},'a*x');
        xFit{windowIdx} = [0,max(meanWindow{windowIdx})*100000];
        if(opts.PLOT_LOG)
            xFit{windowIdx}(1) = eps(1);
        end
        yFit{windowIdx} = F{windowIdx}.a*xFit{windowIdx};
        
        %% get colors if requested
        if(opts.COLOR_MARKERS_RESPONSE && ~isempty(opts.BASELINE_WINDOW_IDX) && windowIdx ~= opts.BASELINE_WINDOW_IDX)
            meanBaseline = meanWindow{opts.BASELINE_WINDOW_IDX};
            if(opts.PROJECT_TO_SLOPE_1) % project data to have a slope of 1 so that mean = variance
                meanBaseline = meanBaseline*F{opts.BASELINE_WINDOW_IDX}.a;
            end
            ratios = meanWindow{windowIdx}./meanBaseline;
            ratios(ratios > opts.MAX_RATIO) = opts.MAX_RATIO;
            ratios(ratios < opts.MIN_RATIO) = opts.MIN_RATIO;
            colorIdx = min(size(opts.COLORS,1),max(1,ceil(size(opts.COLORS,1).*(ratios-opts.MIN_RATIO)./(opts.MAX_RATIO - opts.MIN_RATIO))));
        end
    end
            
    %% plot data for each window
    for windowIdx = 1:numel(meanWindow)
        figHandles{windowIdx} = figure();            
        if(opts.PLOT_LOG && opts.COLOR_MARKERS_RESPONSE && ~isempty(opts.BASELINE_WINDOW_IDX))
            scatter(log10(meanWindow{windowIdx}),log10(varWindow{windowIdx}),opts.MARKER_SIZE,opts.COLORS(colorIdx,:),'filled')
        elseif(opts.PLOT_LOG)
            scatter(log10(meanWindow{windowIdx}),log10(varWindow{windowIdx}),opts.MARKER_SIZE,'b','filled')
        elseif(opts.COLOR_MARKERS_RESPONSE && ~isempty(opts.BASELINE_WINDOW_IDX))
            scatter(meanWindow{windowIdx},varWindow{windowIdx},opts.MARKER_SIZE,opts.COLORS(colorIdx,:),'filled')
        else
            scatter(meanWindow{windowIdx},varWindow{windowIdx},opts.MARKER_SIZE,'b','filled')
        end
        hold on
        
        if(opts.PLOT_FIT_LINE)
            idxUse = windowIdx;
            if(~isempty(opts.BASELINE_WINDOW_IDX))
                idxUse = opts.BASELINE_WINDOW_IDX;
            end
            ax = gca;
            X_LIMITS = ax.XLim;
            if(opts.PLOT_LOG)
                plot(log10(xFit{idxUse}),log10(yFit{idxUse}),'--k','linewidth',2)
            else
                plot(xFit{idxUse},yFit{idxUse},'--k','linewidth',2)
            end
            ax.XLim = X_LIMITS;
        end
        
        globalXLimits = [min(globalXLimits(1),ax.XLim(1)),max(globalXLimits(2),ax.XLim(2))];
        globalYLimits = [min(globalYLimits(1),ax.YLim(1)),max(globalYLimits(2),ax.YLim(2))];
        
        if(opts.PLOT_LOG)
            xlabel('log(Mean)')
            ylabel('log(Variance)')
        else
            xlabel('Mean')
            ylabel('Variance') 
        end

        set(gca,'fontsize',16)
        formatForLee(gcf)
    end
    
    
    %% set limits as the same
    if(opts.SAME_LIMITS)
        for windowIdx = 1:numel(meanWindow)
            figure(figHandles{windowIdx});
            ax = gca;
            ax.XLim = globalXLimits;
            ax.YLim = globalYLimits;
        end
    end
    
    %% deal with saving plots
    for windowIdx = 1:numel(meanWindow)
        if(opts.FIGURE_SAVE && strcmpi(opts.FIGURE_NAME,'')~=1 && strcmpi(opts.FIGURE_DIR,'')~=1)
            strAdd = strcat('_WINDOW_',num2str(opts.WINDOW(windowIdx,1)),'-',num2str(opts.WINDOW(windowIdx,2)));
            if(opts.PROJECT_TO_SLOPE_1)
                strAdd = strcat(strAdd,'_projectedToSlope1');
            end
            if(opts.PLOT_FIT_LINE)
                strAdd = strcat(strAdd,'_fitLine');
            end
            if(opts.PLOT_LOG)
                strAdd = strcat(strAdd,'_logLog');
            end
            saveFiguresLIB(figHandles{windowIdx},opts.FIGURE_DIR,strcat(opts.FIGURE_NAME,strAdd));
        end
    end
    
    % store the fits
    arrayDataFits.F = F;
    arrayDataFits.gof = gof;
    
    %% make colorbar
    figure
    b=colorbar;
    colormap jet;
    set(gca,'Visible','off');
    b.FontSize = 14;
    b.Ticks = [0,0.25,0.5,0.75,1.0];
    maxDataRound = round(opts.MAX_RATIO,1);
    minDataRound = round(opts.MIN_RATIO,1);
    b.TickLabels = {};
    for i = b.Ticks
        if(i==b.Ticks(end))
            b.TickLabels{end+1,1} = strcat('>',num2str(i*(maxDataRound-minDataRound) + minDataRound));
        elseif(i==b.Ticks(1) && i*(maxDataRound-minDataRound) + minDataRound ~= 0)
            b.TickLabels{end+1,1} = strcat('<',num2str(i*(maxDataRound-minDataRound) + minDataRound));
        else
            b.TickLabels{end+1,1} = num2str(i*(maxDataRound-minDataRound) + minDataRound);
        end
    end
    %% deal with saving colorbar plot
    if(opts.FIGURE_SAVE && strcmpi(opts.FIGURE_NAME,'')~=1 && strcmpi(opts.FIGURE_DIR,'')~=1)
        strAdd = strcat('_WINDOW_',num2str(opts.WINDOW(windowIdx,1)),'-',num2str(opts.WINDOW(windowIdx,2)));
        if(opts.PROJECT_TO_SLOPE_1)
            strAdd = strcat(strAdd,'_projectedToSlope1');
        end
        if(opts.PLOT_FIT_LINE)
            strAdd = strcat(strAdd,'_fitLine');
        end
        if(opts.PLOT_LOG)
            strAdd = strcat(strAdd,'_logLog');
        end
        strAdd = strcat(strAdd,'_HEATBAR');
        saveFiguresLIB(gcf,opts.FIGURE_DIR,strcat(opts.FIGURE_NAME,strAdd));
    end
end

function [opts] = configureOpts(optsInput)

    opts.WINDOW = [];
    opts.MARKER_SIZE = 12;
    opts.BASELINE_WINDOW_IDX = [];
    opts.PLOT_FIT_LINE = 1;
    opts.PROJECT_TO_SLOPE_1 = 0;
    opts.COLOR_MARKERS_RESPONSE = 0;
    opts.COLORS = colormap(jet);
    opts.MAX_RATIO = 2;
    opts.MIN_RATIO = 0;
    
    opts.PLOT_LOG = 0;
    opts.SAME_LIMITS = 0;
    
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = '';
    opts.FIGURE_NAME = '';
    
    %% check if in opts and optsInput, overwrite if so
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
    
    if(size(opts.WINDOW,1) == 1)
        opts.BASELINE_WINDOW_IDX = 1;
    end
end