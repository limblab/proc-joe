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
    
    % make figure if requested -- this shows errors as a function of target
    % angle
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
    
    
    %% plot reaches
    % include fail and reward trials, somehow remove outliers that seem
    % like a guess?
    
    % make figure of reaches (all to one fake target or bin based on target?)
    if(opts.MAKE_FIGURES || 1)
        
        % plot all reaches regardless of target direction to start
        timeOffsets = [0,0];
        rotateReaches = 0;
        plotReaches(cds,{'goCueTime','endTime'},timeOffsets,[min(cds.trials.tgtDir),max(cds.trials.tgtDir)],rotateReaches,opts);
        
        
        % plot all reaches to a 0 degree target (combine and rotate)
        timeOffsets = [0,0];
        rotateReaches = 1;
        angleWindow = [-1000,1000];
        plotReaches(cds,{'goCueTime','endTime'},timeOffsets,angleWindow,rotateReaches,opts);
        
        
        % plot distribution of target errors
        angleWindow = [-1000,1000];
        plotReachDistribution(cds,angleWindow,opts);
        
    end
    
    
   
    
end

function [figureHandle] = plotReaches(cds,timeCues,timeOffsets,tgtAngleWindow,rotateReaches,opts)

     % plots reaches between the two time cues offset by time offsets.
     % reaches come from tgts within the tgtAngleWindow. rotateReaches will
     % rotate the space so that all targets are aligned

    if(size(timeOffsets,1) == 1)
        timeOffsets = repmat(timeOffsets,numel(cds.trials.number),1);
    end
     
    figureHandle = figure;
    figureHandle.Position(2) = figureHandle.Position(2)-300;
    figureHandle.Position(4) = figureHandle.Position(3); % make a true circle
    hold on
    numTrialsPlotted = 0;
    
    kinX = cds.kin.x - opts.CENTER_X;
    kinY = cds.kin.y - opts.CENTER_Y; 
    
    for tr = 1:numel(cds.trials.number)
        if((cds.trials.result(tr) == 'R' || cds.trials.result(tr) == 'F') && cds.trials.tgtDir(tr) > tgtAngleWindow(1) && cds.trials.tgtDir(tr) <= tgtAngleWindow(2) && numTrialsPlotted < opts.MAX_TRIALS_PLOT && ...
                ~isnan(cds.trials.(timeCues{1})(tr)) && ~isnan(cds.trials.(timeCues{2})(tr)))
            kinIdx = [find(cds.trials.(timeCues{1})(tr)+timeOffsets(tr,1) <= cds.kin.t,1,'first'),find(cds.trials.(timeCues{2})(tr)+timeOffsets(tr,2) <= cds.kin.t,1,'first')];
            xData = kinX(kinIdx(1):kinIdx(2));
            yData = kinY(kinIdx(1):kinIdx(2));

            if(rotateReaches)
                xDataRot = xData*cos(cds.trials.tgtDir(tr)*pi/180) + yData*sin(cds.trials.tgtDir(tr)*pi/180);
                yDataRot = -xData*sin(cds.trials.tgtDir(tr)*pi/180) + yData*cos(cds.trials.tgtDir(tr)*pi/180);
                xData = xDataRot;
                yData = yDataRot;
            end
            if(opts.PLOT_UNTIL_CIRCLE)
                outerCircleIdx = find(sqrt((xData.^2 + yData.^2)) > opts.CIRCLE_RADIUS - 0.5*opts.CIRCLE_DEPTH,1,'first');
                xData = xData(1:outerCircleIdx-1);
                yData = yData(1:outerCircleIdx-1);
            end
            plot(xData,yData)
            numTrialsPlotted = numTrialsPlotted + 1;
        end
    end
     
    if(opts.DRAW_CIRCLE)
        t = linspace(0,2*pi,10000);
        x = (opts.CIRCLE_RADIUS-0.5*opts.CIRCLE_DEPTH)*cos(t);
        y = (opts.CIRCLE_RADIUS-0.5*opts.CIRCLE_DEPTH)*sin(t);
        plot(x,y,'k','linewidth',2)
    end
    
    xlim(opts.XLIM)
    ylim(opts.YLIM)
    
    if(opts.REMOVE_AXIS)
        ax = gca;
        set(ax,'Visible','off');
    end
    
end


function [figureHandle] = plotReachDistribution(cds,tgtAngleWindow,opts)

    % plots distribution of final reach angle relative to target direction
    figureHandle = figure;
    
    kinX = cds.kin.x - opts.CENTER_X;
    kinY = cds.kin.y - opts.CENTER_Y;
    
    reachAngles = zeros(numel(cds.trials.number),1);
    sizeReachAngle = 0;
    angleOffsets = [-360,0,360];
    for tr = 1:numel(cds.trials.number)
        if((cds.trials.result(tr) == 'R' || cds.trials.result(tr) == 'F') && cds.trials.tgtDir(tr) > tgtAngleWindow(1) && cds.trials.tgtDir(tr) <= tgtAngleWindow(2) && ...
                ~isnan(cds.trials.goCueTime(tr)) && ~isnan(cds.trials.endTime(tr)))
            kinIdx = [find(cds.trials.goCueTime(tr) <= cds.kin.t,1,'first'),find(cds.trials.endTime(tr) <= cds.kin.t,1,'first')];
            xData = kinX(kinIdx(1):kinIdx(2));
            yData = kinY(kinIdx(1):kinIdx(2));
            
            outerCircleIdx = find(sqrt((xData.^2 + yData.^2)) > opts.CIRCLE_RADIUS - 0.5*opts.CIRCLE_DEPTH - 0.5,1,'first');
            reachAngle = 180/pi*atan2(yData(outerCircleIdx-1),xData(outerCircleIdx-1));
            [~,offsetIdx] = min(abs(cds.trials.tgtDir(tr) - (reachAngle+angleOffsets)));
            reachAngles(sizeReachAngle+1,1) = cds.trials.tgtDir(tr) - (reachAngle+angleOffsets(offsetIdx));
            sizeReachAngle = sizeReachAngle + 1;
        end
    end

    reachAngles = reachAngles(1:sizeReachAngle);
    % bin and plot reach angle
    
    [bC,bE] = histcounts(reachAngles,'BinWidth',opts.DISTRIBUTION_BIN_SIZE);
    bC = bC/numel(reachAngles);
    bar(bE(1:end-1) + (bE(2)-bE(1))/2,bC);
    
    formatForLee(gcf)
end


function [opts] = configureOpts(optsInput)

    opts = [];
    opts.NUM_BINS_DIR = 8;
    opts.MAKE_FIGURES = 1;
    opts.MARKER_SIZE = 12;
    opts.PLOT_POLAR = 0;
    opts.PLOT_CHANCE = 1;
    opts.FONT_SIZE = 16;
    
    opts.CENTER_X = 3;
    opts.CENTER_Y = -33;
    
    opts.MAX_TRIALS_PLOT = 100;
    opts.DRAW_CIRCLE = 1;
    opts.CIRCLE_RADIUS = 7;
    opts.CIRCLE_DEPTH = 2;
    opts.PLOT_UNTIL_CIRCLE = 1;
    
    opts.XLIM = [-opts.CIRCLE_RADIUS, opts.CIRCLE_RADIUS];
    opts.YLIM = [-opts.CIRCLE_RADIUS, opts.CIRCLE_RADIUS];
    
    opts.REMOVE_AXIS = 1;
    
    opts.DISTRIBUTION_BIN_SIZE = 5;
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