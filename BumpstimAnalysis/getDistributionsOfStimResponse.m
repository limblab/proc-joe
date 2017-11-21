function [arrayDataFits] = getDistributionsOfStimResponse( arrayData,opts )
% 

    opts = configureOpts(opts);

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
    end % end for window
    
    %% plot data for each window with a best fit line going through (0,0)
    windowIdxs = 1:1:size(opts.WINDOW,1);
    if(~isempty(opts.BASELINE_WINDOW_IDX))
        windowIdxs(opts.BASELINE_WINDOW_IDX) = 1;
        windowIdxs(1) = opts.BASELINE_WINDOW_IDX;
    end
    for windowIdx = windowIdxs
        % get data and compress across all units, channels, waveforms sent
        meanWindow = reshape(means(windowIdx,:,:,:),numel(means(windowIdx,:,:,:)),1);
        varWindow = reshape(vars(windowIdx,:,:,:),numel(vars(windowIdx,:,:,:)),1);
        
        % plot variance vs mean
        
        [F{windowIdx},gof{windowIdx}] = fit(meanWindow(meanWindow < 0.75*max(meanWindow)),varWindow(meanWindow < 0.75*max(meanWindow)),'a*x');
        xFit = [0,max(meanWindow)];
        if(opts.PLOT_LOG)
            xFit(1) = eps(1);
        end
        
        if(isempty(opts.BASELINE_WINDOW_IDX))
            yFit = F{windowIdx}.a*xFit;
        else
            yFit = F{opts.BASELINE_WINDOW_IDX}.a*xFit;
        end
        
        if(~isempty(opts.BASELINE_WINDOW_IDX) && opts.PROJECT_TO_SLOPE_1)
            xFit = xFit*F{opts.BASELINE_WINDOW_IDX}.a;
            meanWindow = meanWindow*F{opts.BASELINE_WINDOW_IDX}.a;
        end
        figure();
        if(opts.COLOR_MARKERS_RESPONSE && ~isempty(opts.BASELINE_WINDOW_IDX) && windowIdx ~= opts.BASELINE_WINDOW_IDX)
            meanBaseline = reshape(means(opts.BASELINE_WINDOW_IDX,:,:,:),numel(means(opts.BASELINE_WINDOW_IDX,:,:,:)),1);
            if(opts.PROJECT_TO_SLOPE_1)
                meanBaseline = meanBaseline*F{opts.BASELINE_WINDOW_IDX}.a;
            end
            ratios = meanWindow./meanBaseline;
            ratios(ratios > opts.MAX_RATIO) = opts.MAX_RATIO;
            ratios(ratios < opts.MIN_RATIO) = opts.MIN_RATIO;
            
            for r = 1:numel(ratios)
                colorIdx = min(size(opts.COLORS,1),max(1,ceil(size(opts.COLORS,1)*(ratios(r)-opts.MIN_RATIO)/(opts.MAX_RATIO - opts.MIN_RATIO))));
                
                if(opts.PLOT_LOG)
                    plot(log(meanWindow(r)),log(varWindow(r)),'.','markersize',opts.MARKER_SIZE,'color',opts.COLORS(colorIdx,:));
                else
                    plot(meanWindow(r),varWindow(r),'.','markersize',opts.MARKER_SIZE,'color',opts.COLORS(colorIdx,:));
                end
                hold on
            end
        else
            if(opts.PLOT_LOG)
                plot(log(meanWindow),log(varWindow),'.','markersize',opts.MARKER_SIZE);
            else
                plot(meanWindow,varWindow,'.','markersize',opts.MARKER_SIZE);
            end
        end
        hold on
        if(opts.PLOT_FIT_LINE)
            ax = gca;
            X_LIMITS = ax.XLim;
            if(opts.PLOT_LOG)
                plot(log(xFit),log(yFit),'--k','linewidth',2)
            else
                plot(xFit,yFit,'--k','linewidth',2)
            end
            ax.XLim = X_LIMITS;
        end
        if(opts.PLOT_LOG)
            xlabel('log(Mean)')
            ylabel('log(Variance)')
        else
            xlabel('Mean')
            ylabel('Variance') 
        end

        set(gca,'fontsize',16)
        formatForLee(gcf)
        
        %% deal with saving plot
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
            saveFiguresLIB(gcf,opts.FIGURE_DIR,strcat(opts.FIGURE_NAME,strAdd));
        end
        
    end
    
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
    %% deal with saving plot
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