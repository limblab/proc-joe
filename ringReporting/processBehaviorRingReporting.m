function [behaviorData,figureHandles] = processBehaviorRingReporting(cds,opts)

    
    %% configure opts
    opts = configureOpts(opts);
    figureHandles = {};
    
    %% store % correct vs % fail
    if(isempty(opts.BUMP_MAGS))
        opts.BUMP_MAGS = unique(cds.trials.bumpMagnitude);
        opts.BUMP_MAGS = opts.BUMP_MAGS(~isnan(opts.BUMP_MAGS));
    end
    
    for bm = 1:numel(opts.BUMP_MAGS)
        behaviorData.percentCorrect(bm) = sum(cds.trials.result == 'R' & cds.trials.hideCursorDuringBump == 1 & cds.trials.showOuterTarget == 0 & ~isnan(cds.trials.bumpMagnitude) & cds.trials.bumpMagnitude == opts.BUMP_MAGS(bm))/...
            sum((cds.trials.result == 'R' | cds.trials.result == 'F') & cds.trials.hideCursorDuringBump == 1 & cds.trials.showOuterTarget == 0 & ~isnan(cds.trials.bumpMagnitude) & cds.trials.bumpMagnitude == opts.BUMP_MAGS(bm));
    end
    
    behaviorData.percentCorrect_all = sum(cds.trials.result == 'R' & cds.trials.hideCursorDuringBump == 1 & cds.trials.showOuterTarget == 0 & ~isnan(cds.trials.bumpMagnitude))/...
            sum((cds.trials.result == 'R' | cds.trials.result == 'F') & cds.trials.hideCursorDuringBump == 1 & cds.trials.showOuterTarget == 0 & ~isnan(cds.trials.bumpMagnitude));
    
%     %% look at % correct as a function of trials
%     trials_attempted = find((cds.trials.result == 'R' | cds.trials.result == 'F') & cds.trials.showOuterTarget == 0 & ~isnan(cds.trials.bumpMagnitude));
%     trials_result = cds.trials.result(trials_attempted) == 'R';
%     
%     for i = 1:floor(numel(trials_result)/20)
%         binned_result(i) = sum(trials_result((i-1)*20+1:(i-1)*20+20));
%     end
    %% look at % correct as a function of tgtDir (make plot) and bump mag
    % opts.NUM_BINS_DIR is the number of bins to divide the circle into
     
    behaviorData.bin_edges = (0:360/opts.NUM_BINS_DIR:360) - 180;
    
    behaviorData.bin_percentCorrect = zeros(1,opts.NUM_BINS_DIR);
    for bm = 1:numel(opts.BUMP_MAGS)
        for b = 1:opts.NUM_BINS_DIR
            behaviorData.bin_percentCorrect(bm,b) = sum(cds.trials.result == 'R' & cds.trials.hideCursorDuringBump == 1 & cds.trials.showOuterTarget == 0 & cds.trials.bumpMagnitude == opts.BUMP_MAGS(bm) &...
                cds.trials.tgtDir >= behaviorData.bin_edges(b) & cds.trials.tgtDir < behaviorData.bin_edges(b+1))/...
                sum((cds.trials.result == 'R' | cds.trials.result == 'F') & cds.trials.hideCursorDuringBump == 1 & cds.trials.showOuterTarget == 0 & cds.trials.bumpMagnitude == opts.BUMP_MAGS(bm) &...
                cds.trials.tgtDir >= behaviorData.bin_edges(b) & cds.trials.tgtDir < behaviorData.bin_edges(b+1));
        end
    end
    
    % make figure if requested -- this shows errors as a function of target
    % angle
    if(opts.MAKE_FIGURES && opts.NUM_BINS_DIR > 1)
        figureHandles{end+1} = figure;
        for bm = 1:numel(opts.BUMP_MAGS)
            if(opts.PLOT_POLAR)
                theta = deg2rad([behaviorData.bin_edges(1:end-1),behaviorData.bin_edges(1)] + (behaviorData.bin_edges(2)-behaviorData.bin_edges(1))/2);
                rho = [behaviorData.bin_percentCorrect(bm,:),behaviorData.bin_percentCorrect(bm,1)]; % closing the cirlce by including the first value twice
                polarplot(theta,rho,'-.','markersize',opts.MARKER_SIZE,'linewidth',opts.LINE_WIDTH,'color',opts.COLORS(bm,:));
                figureHandles{end}.Children(1).RLim = [0,1]; % radius limits from 0 to 1
                ax = gca;
                ax.ThetaTickLabel = {'0','30','60','90','120','150','180','-150','-120','-90','-60','-30'};
            else % line plot
                plot(behaviorData.bin_edges(1:end-1) + (behaviorData.bin_edges(2)-behaviorData.bin_edges(1))/2,behaviorData.bin_percentCorrect(bm,:),'-.','markerSize',opts.MARKER_SIZE,'linewidth',opts.LINE_WIDTH);
                ylim([0,1]);
                xlim([-180,180]);
                set(gca,'fontsize',opts.FONT_SIZE)
                formatForLee(gcf)
            end
            set(gca,'fontsize',14);
            hold on
        end
        
        if(opts.FIGURE_SAVE)
            saveFiguresLIB(gcf,opts.FIGURE_DIR,strcat(opts.FIGURE_PREFIX,'_percentCorrect'));
        end
    end
    
    
    %% plot reaches
    % include fail and reward trials, somehow remove outliers that seem
    % like a guess?
    
    % make figure of reaches (all to one fake target or bin based on target?)
    if(opts.MAKE_FIGURES || 1)

%         % plot all reaches regardless of target direction to start
%         timeOffsets = [0,0];
%         rotateReaches = 0;
%         removeBumpOffset = 0;
%         [~,reachData_withOffset] = plotReaches(cds,{'goCueTime','endTime'},timeOffsets,[min(cds.trials.tgtDir),max(cds.trials.tgtDir)],rotateReaches,removeBumpOffset,opts);
%         
        
%         for i = 1:numel(reachData_normal.x)
%             try
%                 if(~(atan2(reachData_normal.y{i}(end),reachData_normal.x{i}(end)) < pi && ...
%                         atan2(reachData_normal.y{i}(end),reachData_normal.x{i}(end)) > pi/2))
%                     diff(i) = acosd(sum([reachData_normal.x{i}(end) reachData_normal.y{i}(end)].*...
%                         [reachData_withOffset.x{i}(end) reachData_withOffset.y{i}(end)])./...
%                         (sqrt(sum([reachData_normal.x{i}(end) reachData_normal.y{i}(end)].^2))*...
%                         sqrt(sum([reachData_withOffset.x{i}(end) reachData_withOffset.y{i}(end)].^2))));
%                 else
%                     diff(i) = -1;
%                 end
%             catch
%             end
%         end
        
% %         % plot all reaches to a 0 degree target (combine and rotate)
%         timeOffsets = [0,0];
%         rotateReaches = 1;
%         angleWindow = [-1000,1000];
%         removeBumpOffset = 0;
%         plotReaches(cds,{'goCueTime','endTime'},timeOffsets,angleWindow,rotateReaches,removeBumpOffset,opts);

        % plot distribution of target errors
        angleWindow = [-1000,1000];
        opts.PLOT_ALL_BUMP_MAGS = 1;
        removeBumpOffset = 0;
        plotReachDistribution(cds,angleWindow,removeBumpOffset,opts);

%         % plot kinematics during a bump
%         timeOffsets = [0,0];
%         angleWindow = [-150,-50];
%         rotateReaches = 0;
%         removeBumpOffset = 1;
%         [~,~] = plotReaches(cds,{'bumpTime','goCueTime'},timeOffsets,angleWindow,rotateReaches,removeBumpOffset,opts);
%         [~,reaches] = plotReaches(cds,{'goCueTime','endTime'},timeOffsets,angleWindow,rotateReaches,removeBumpOffset,opts);
%         reach_angle = zeros(numel(reaches),1);
%         for ra = 1:numel(reaches.x)
%             reach_angle(ra) = 180/pi*atan2(reaches.y{ra}(end)-reaches.y{ra}(1),reaches.x{ra}(end)-reaches.x{ra}(1));
%         end
%         figure
%         [bC,bE] = histcounts(reach_angle,'BinWidth',opts.DISTRIBUTION_BIN_SIZE);
%         bar(bE(1:end-1) + mode(diff(bE))/2,bC);
    end
    
%     %% plot polar plot showing distribution of target centers to check for randomness
%     [tgtDistribution] = histcounts(cds.trials.tgtDir((cds.trials.result == 'R' | cds.trials.result == 'F') & ~isnan(cds.trials.tgtDir)),behaviorData.bin_edges);
%     figureHandles{end+1} = figure;
%     theta = deg2rad([behaviorData.bin_edges(1:end-1),behaviorData.bin_edges(1)] + (behaviorData.bin_edges(2)-behaviorData.bin_edges(1))/2);
%     rho = [tgtDistribution,tgtDistribution(1)]/sum((cds.trials.result == 'R' | cds.trials.result == 'F') & ~isnan(cds.trials.tgtDir)); % closing the cirlce by including the first value twice
%     polarplot(theta,rho,'-.','markersize',opts.MARKER_SIZE,'linewidth',opts.LINE_WIDTH);
%     title('target distribution for completed reaches')
%     figureHandles{end}.Children(1).RLim = [0,1/opts.NUM_BINS_DIR + 0.1]; % radius limits from 0 to 1
%     
%     
   
    
end

function [figureHandle,reachData] = plotReaches(cds,timeCues,timeOffsets,tgtAngleWindow,rotateReaches,removeBumpOffset,opts)

     % plots reaches between the two time cues offset by time offsets.
     % reaches come from tgts within the tgtAngleWindow. rotateReaches will
     % rotate the space so that all targets are aligned   
     reachData.x = {};
     reachData.y = {};
     
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
                ~isnan(cds.trials.(timeCues{1})(tr)) && ~isnan(cds.trials.(timeCues{2})(tr)) && cds.trials.hideCursorDuringBump(tr) == 1 && cds.trials.showRing(tr) && ~cds.trials.showOuterTarget(tr))
            kinIdx = [find(cds.trials.(timeCues{1})(tr)+timeOffsets(tr,1) <= cds.kin.t,1,'first'),find(cds.trials.(timeCues{2})(tr)+timeOffsets(tr,2) <= cds.kin.t,1,'first')];
            xData = kinX(kinIdx(1):kinIdx(2));
            yData = kinY(kinIdx(1):kinIdx(2));
            
            if(removeBumpOffset)
                xData = xData - kinX(kinIdx(1));
                yData = yData - kinY(kinIdx(1));
            end
            if(rotateReaches)
                xDataRot = xData*cos(cds.trials.tgtDir(tr)*pi/180) + yData*sin(cds.trials.tgtDir(tr)*pi/180);
                yDataRot = -xData*sin(cds.trials.tgtDir(tr)*pi/180) + yData*cos(cds.trials.tgtDir(tr)*pi/180);
                xData = xDataRot;
                yData = yDataRot;
            end
            if(opts.PLOT_UNTIL_CIRCLE)
                outerCircleIdx = find(sqrt((xData.^2 + yData.^2)) > opts.CIRCLE_RADIUS - 0.5*opts.CIRCLE_DEPTH,1,'first');
                if(~isempty(outerCircleIdx))
                    xData = xData(1:outerCircleIdx-1);
                    yData = yData(1:outerCircleIdx-1);
                end
            end
            plot(xData,yData)
            reachData.x{end+1} = xData;
            reachData.y{end+1} = yData;
            numTrialsPlotted = numTrialsPlotted + 1;
        end
    end
     
    if(opts.DRAW_CIRCLE)
        t = linspace(0,2*pi,10000);
        x = (opts.CIRCLE_RADIUS-0.5*opts.CIRCLE_DEPTH)*cos(t);
        y = (opts.CIRCLE_RADIUS-0.5*opts.CIRCLE_DEPTH)*sin(t);
        plot(x,y,'k','linewidth',opts.CIRCLE_LINE_WIDTH)
        hold on
        if(rotateReaches)
            t = linspace(-cds.trials.tgtWidth(1)/2,cds.trials.tgtWidth(1)/2,10000)*pi/180;
            x = (opts.CIRCLE_RADIUS-0.5*opts.CIRCLE_DEPTH)*cos(t);
            y = (opts.CIRCLE_RADIUS-0.5*opts.CIRCLE_DEPTH)*sin(t);
            plot(x,y,'r','linewidth',opts.CIRCLE_LINE_WIDTH)
        else
            t = linspace(-cds.trials.tgtWidth(1)/2,cds.trials.tgtWidth(1)/2,10000)*pi/180;
            t = t + mean(tgtAngleWindow)*pi/180;
            x = (opts.CIRCLE_RADIUS-0.5*opts.CIRCLE_DEPTH)*cos(t);
            y = (opts.CIRCLE_RADIUS-0.5*opts.CIRCLE_DEPTH)*sin(t);
            plot(x,y,'r','linewidth',opts.CIRCLE_LINE_WIDTH)
        end
        
        % plot center circle
        plot(0,0,'k.','markersize',10)
    end
    
    xlim(opts.XLIM)
    ylim(opts.YLIM)
    
    if(opts.REMOVE_AXIS)
        ax = gca;
        set(ax,'Visible','off');
    end
    
    if(opts.FIGURE_SAVE)
        if(rotateReaches)
            saveFiguresLIB(gcf,opts.FIGURE_DIR,strcat(opts.FIGURE_PREFIX,'_reachesRotate'));
        else
            saveFiguresLIB(gcf,opts.FIGURE_DIR,strcat(opts.FIGURE_PREFIX,'_reaches'));
        end
    end
end

function [figureHandle] = plotReachDistribution(cds,tgtAngleWindow,removeBumpOffset,opts)

    % plots distribution of final reach angle relative to target direction
    figureHandle = figure;
    
    kinX = cds.kin.x - opts.CENTER_X;
    kinY = cds.kin.y - opts.CENTER_Y;
    
    reachAngles = zeros(numel(cds.trials.number),1);
    reachError = zeros(numel(cds.trials.number),1);
    bumpAngles = zeros(numel(cds.trials.number),1);
    bumpMag = zeros(numel(cds.trials.number),1);
    sizeReachAngle = 0;
    angleOffsets = [-360,0,360];
    
    
    for tr = 1:numel(cds.trials.number)
        if((cds.trials.result(tr) == 'R' || cds.trials.result(tr) == 'F') && cds.trials.tgtDir(tr) > tgtAngleWindow(1) && cds.trials.tgtDir(tr) <= tgtAngleWindow(2) && ...
                ~isnan(cds.trials.goCueTime(tr)) && ~isnan(cds.trials.endTime(tr)) && cds.trials.hideCursorDuringBump(tr) == 1 && cds.trials.showRing(tr) && ~cds.trials.showOuterTarget(tr))
            kinIdx = [find(cds.trials.goCueTime(tr) <= cds.kin.t,1,'first'),find(cds.trials.endTime(tr) <= cds.kin.t,1,'first')];
            xData = kinX(kinIdx(1):kinIdx(2));
            yData = kinY(kinIdx(1):kinIdx(2));
            tData = cds.kin.t(kinIdx(1):kinIdx(2));
            
            if(removeBumpOffset)
                xData = xData - kinX(kinIdx(1));
                yData = yData - kinY(kinIdx(1));
            end
            
            if(ismember('otHoldTime',cds.trials.Properties.VariableNames)) % 
                outerCircleIdx = find(tData > cds.trials.otHoldTime(tr),1,'first');
            else
                outerCircleIdx = find(sqrt((xData.^2 + yData.^2)) > opts.CIRCLE_RADIUS - 0.5*opts.CIRCLE_DEPTH,1,'first');
                if(isempty(outerCircleIdx)) 
                    outerCircleIdx = numel(xData);
                end
            end
            reachAngles(sizeReachAngle+1,1) = 180/pi*atan2(yData(outerCircleIdx-1),xData(outerCircleIdx-1));
            [~,offsetIdx] = min(abs(cds.trials.tgtDir(tr) - (reachAngles(sizeReachAngle+1,1)+angleOffsets)));
            reachError(sizeReachAngle+1,1) = cds.trials.tgtDir(tr) - (reachAngles(sizeReachAngle+1,1)+angleOffsets(offsetIdx));
            bumpAngles(sizeReachAngle+1,1) = cds.trials.bumpDir(tr);
            bumpMag(sizeReachAngle+1,1) = cds.trials.bumpMagnitude(tr);
            sizeReachAngle = sizeReachAngle + 1;
        end
    end

    reachAngles = reachAngles(1:sizeReachAngle,1);
    reachError = reachError(1:sizeReachAngle);
    bumpAngles = bumpAngles(1:sizeReachAngle);
    bumpMag = bumpMag(1:sizeReachAngle);
    
    if(opts.PLOT_ALL_BUMP_MAGS)
        for bm = 1:numel(opts.BUMP_MAGS)
            [bC,bE] = histcounts(reachError(bumpMag == opts.BUMP_MAGS(bm)),'BinWidth',opts.DISTRIBUTION_BIN_SIZE);
            bC = bC/sum(bumpMag == opts.BUMP_MAGS(bm));
            plot(bE(1:end-1) + (bE(2)-bE(1))/2,bC,'linewidth',opts.LINE_WIDTH,'color',opts.COLORS(bm,:));
            hold on
        end
    else
        % bin and plot reach angle
        [bC,bE] = histcounts(reachError,'BinWidth',opts.DISTRIBUTION_BIN_SIZE);
        bC = bC/numel(reachError);
        bar(bE(1:end-1) + (bE(2)-bE(1))/2,bC);
        xlabel('Error to target (deg)')
        ylabel('Proportion of reaches')
        ax = gca;
        set(ax,'fontsize',14)
        hold on
        plot(-1*mode(cds.trials.tgtWidth)*[1/2,1/2],ax.YLim,'r--','linewidth',2)
        plot(mode(cds.trials.tgtWidth)*[1/2,1/2],ax.YLim,'r--','linewidth',2)
    end
    formatForLee(gcf)
    
    if(opts.FIGURE_SAVE)
        saveFiguresLIB(gcf,opts.FIGURE_DIR,strcat(opts.FIGURE_PREFIX,'_errorToTarget'));
    end
    
    % plot reach direction vs bump direction
    figure();
    if(opts.PLOT_ALL_BUMP_MAGS)
        for bm = 1:numel(opts.BUMP_MAGS)
            plot(bumpAngles(abs(bumpAngles) < 180 & abs(reachAngles) < 180 & bumpMag == opts.BUMP_MAGS(bm)),...
                reachAngles(abs(bumpAngles) < 180 & abs(reachAngles) < 180 & bumpMag == opts.BUMP_MAGS(bm)),...
                '.','markersize',10,'color',opts.COLORS(bm,:));
            hold on
        end
    else
        plot(bumpAngles(abs(bumpAngles) < 180 & abs(reachAngles) < 180),...
            reachAngles(abs(bumpAngles) < 180 & abs(reachAngles) < 180),...
            '.','markersize',10);
    end
    xlabel('Bump direction (deg)')
    ylabel('Reach direction (deg)')
    hold on
    xDeg = -180:1:180;
    plot(xDeg,xDeg,'k--','linewidth',2) % plots unity
    plot(xDeg,xDeg + mode(cds.trials.tgtWidth)/2,'r--','linewidth',2);
    plot(xDeg,xDeg - mode(cds.trials.tgtWidth)/2,'r--','linewidth',2);
    % deal with wrap around
    plot(xDeg,xDeg + mode(cds.trials.tgtWidth)/2-360,'r--','linewidth',2);
    plot(xDeg,xDeg - mode(cds.trials.tgtWidth)/2+360,'r--','linewidth',2);
    ylim([-180,180])
    xlim([-180,180])
    set(gca,'fontsize',14);
    formatForLee(gcf)
    
    if(opts.FIGURE_SAVE)
        saveFiguresLIB(gcf,opts.FIGURE_DIR,strcat(opts.FIGURE_PREFIX,'_bumpVsReach'));
    end
    
end


function [opts] = configureOpts(optsInput)

    opts = [];
    opts.NUM_BINS_DIR = 8;
    opts.MAKE_FIGURES = 1;
    opts.MARKER_SIZE = 12;
    opts.PLOT_POLAR = 0;
    opts.PLOT_CHANCE = 1;
    opts.FONT_SIZE = 14;
    opts.LINE_WIDTH = 2;
    opts.CENTER_X = 3;
    opts.CENTER_Y = -33;
    
    opts.MAX_TRIALS_PLOT = 100;
    opts.DRAW_CIRCLE = 1;
    opts.CIRCLE_RADIUS = 7;
    opts.CIRCLE_DEPTH = 2;
    opts.PLOT_UNTIL_CIRCLE = 1;
    opts.CIRCLE_LINE_WIDTH = 3;
    
    opts.XLIM = [-opts.CIRCLE_RADIUS, opts.CIRCLE_RADIUS];
    opts.YLIM = [-opts.CIRCLE_RADIUS, opts.CIRCLE_RADIUS];
    
    opts.REMOVE_AXIS = 1;
    opts.BUMP_MAGS = [];
    
    opts.DISTRIBUTION_BIN_SIZE = 5;
    
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_PREFIX = '';
    opts.FIGURE_DIR = '';
    
    opts.COLORS = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]/255;
    
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