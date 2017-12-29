function [behaviorData,figureHandles] = processBehaviorRingReporting(cds,opts)

    
    %% configure opts
    opts = configureOpts(opts);
    figureHandles = {};
    
    %% store % correct vs % fail
    behaviorData.percentCorrect = sum(cds.trials.result == 'R')/sum(cds.trials.result == 'R' | cds.trials.result == 'F');
    
    %% look at % correct as a function of tgtDir (make plot)
    % opts.NUM_BINS_DIR is the number of bins to divide the circle into
    
    behaviorData.bin_edges = (0:360/opts.NUM_BINS_DIR:360-360/opts.NUM_BINS_DIR) - (360/opts.NUM_BINS_DIR)/2 - 180;
    behaviorData.bin_edges(end+1) = behaviorData.bin_edges(1) + 360;
    
    
    
    behaviorData.bin_percentCorrect = zeros(1,opts.NUM_BINS_DIR);
    for b = 1:opts.NUM_BINS_DIR
        behaviorData.bin_percentCorrect(b) = sum(cds.trials.result == 'R' & ...
            cds.trials.tgtDir >= behaviorData.bin_edges(b) & cds.trials.tgtDir < behaviorData.bin_edges(b+1))/...
            sum((cds.trials.result == 'R' | cds.trials.result == 'F') & cds.trials.tgtDir >= behaviorData.bin_edges(b) & cds.trials.tgtDir < behaviorData.bin_edges(b+1));
    end
    
    % make figure if requested
    if(opts.MAKE_FIGURES && opts.NUM_BINS_DIR > 1)
        figureHandles{end+1} = figure;
        if(opts.PLOT_POLAR)
            theta = deg2rad([behaviorData.bin_edges(1:end-1),behaviorData.bin_edges(1)] + (behaviorData.bin_edges(2)-behaviorData.bin_edges(1))/2);
            rho = [behaviorData.bin_percentCorrect,behaviorData.bin_percentCorrect(1)]; % closing the cirlce by including the first value twice
            polarplot(theta,rho,'-k.','markersize',opts.MARKER_SIZE);
            figureHandles{end}.Children(1).RLim = [0,1]; % radius limits from 0 to 1
        else % line plot
            plot(behaviorData.bin_edges(1:end-1) + (behaviorData.bin_edges(2)-behaviorData.bin_edges(1))/2,behaviorData.bin_percentCorrect,'-k.','markerSize',opts.MARKER_SIZE);
            ylim([0,1]);
            xlim([-180,180]);
            set(gca,'fontsize',opts.FONT_SIZE)
            formatForLee(gcf)
        end
        
    end
    
    
    %% get an error metric from the center of a target
    % include fail and reward trials, somehow remove outliers that seem
    % like a guess?
    
    % make figure of reaches (all to one fake target or bin based on target?)
    if(opts.MAKE_FIGURES || 1)
        % plot all reaches as if target was at 0 degrees
        figureHandles{end+1} = figure;
        hold on
        for tr = 1:numel(cds.trials.number)
            if(cds.trials.result(tr) == 'R')
                kinIdx = [find(round(cds.trials.goCueTime(tr)-0.3,2) == cds.kin.t),find(round(cds.trials.endTime(tr)-0.3,2) == cds.kin.t)];
                xData = cds.kin.x(kinIdx(1):kinIdx(2));
                yData = cds.kin.y(kinIdx(1):kinIdx(2));
                plot(xData,yData)
            end
        end
    end
end


function [opts] = configureOpts(optsInput)

    opts = [];
    opts.NUM_BINS_DIR = 8;
    opts.MAKE_FIGURES = 1;
    opts.MARKER_SIZE = 12;
    opts.PLOT_POLAR = 0;
    opts.PLOT_CHANCE = 1;
    opts.FONT_SIZE = 16;
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
end