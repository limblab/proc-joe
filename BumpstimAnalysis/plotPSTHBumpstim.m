function [ figureHandle ] = plotPSTHBumpstim( cds,neuronNumber,optsTaskInput,optsPlotInput,optsSaveInput )
% this function plots rasters for the bumpstim task(s) based around the
% trial data
%% configure options and initialize output
    optsPlot = configureOptionsPlot(optsPlotInput);
    optsSave = configureOptionsSave(optsSaveInput);
    optsTask = configureOptionsTask(optsTaskInput,cds);
    
    figureHandle = '';
    
    
    %% set nan's in cds.trials.stimCode to -1 and bumpDir
    cds.trials.stimCode(isnan(cds.trials.stimCode)) = -1;
    cds.trials.bumpDir(isnan(cds.trials.bumpDir)) = -1;
    %% set tgtDir, bumpDir, stimCode to -100 if in optsTask.IGNORE
    for ig = 1:numel(optsTask.IGNORE)
        switch optsTask.IGNORE{ig}
            case 'tgtDir'
                optsTask.TARGET_DIRECTIONS = -100;
                cds.trials.tgtDir(~isnan(cds.trials.tgtDir)) = -100;
            case 'bumpDir'
                optsTask.BUMP_DIRECTIONS = -100;
                cds.trials.bumpDir(~isnan(cds.trials.bumpDir)) = -100;
        end
    end
    %% extract data
    spikeTimeData = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),numel(optsTask.BUMP_DIRECTIONS),numel(optsTask.STIM_CODE));
    trialCounter = zeros(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),numel(optsTask.BUMP_DIRECTIONS),numel(optsTask.STIM_CODE));
    if(isfield(cds,'stimOn'))
        stimData = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),numel(optsTask.BUMP_DIRECTIONS),numel(optsTask.STIM_CODE));
    end
    
    stimIdxStart = 1;
    stimIdxEnd = 100;
    binMin = -1;
    binMax = 1;
    
    for trialNum = 1:size(cds.trials,1)
        if(cds.trials.result(trialNum) == 'R' && ~isnan(cds.trials.tgtOnTime(trialNum)) && ~isnan(cds.trials.movePeriod(trialNum))) % if a reward trial and not a corrupted trial
            taskString = '';
            
            if(cds.trials.ctrHoldBump(trialNum) && cds.trials.stimCode(trialNum)==-1 && sum(strcmpi(optsTask.TRIAL_LIST,'ctrHoldBump')) > 0)
                taskString = 'ctrHoldBump';
            elseif(cds.trials.ctrHoldBump(trialNum) && cds.trials.stimCode(trialNum)~=-1 && sum(strcmpi(optsTask.TRIAL_LIST,'ctrHoldBumpStim')) > 0)
                taskString = 'ctrHoldBump';
            elseif(cds.trials.delayBump(trialNum) && cds.trials.stimCode(trialNum)==-1 && sum(strcmpi(optsTask.TRIAL_LIST,'delayBump')) > 0)
                taskString = 'delayBump';
            elseif(cds.trials.delayBump(trialNum) && cds.trials.stimCode(trialNum)~=-1 && sum(strcmpi(optsTask.TRIAL_LIST,'delayBumpStim')) > 0)
                taskString = 'delayBump';
            elseif(cds.trials.moveBump(trialNum) && cds.trials.stimCode(trialNum)==-1 && sum(strcmpi(optsTask.TRIAL_LIST,'moveBump')) > 0)
                taskString = 'moveBump';
            elseif(cds.trials.moveBump(trialNum) && cds.trials.stimCode(trialNum)~=-1 && sum(strcmpi(optsTask.TRIAL_LIST,'moveBumpStim')) > 0)
                taskString = 'moveBump';
            elseif(~cds.trials.ctrHoldBump(trialNum) && ~cds.trials.delayBump(trialNum) && ~cds.trials.moveBump(trialNum) && sum(strcmpi(optsTask.TRIAL_LIST,'noBump')) > 0)
                taskString = 'noBump';
            end
            
            if(strcmpi(taskString,'')~=1 && (~isempty(find(optsTask.TARGET_DIRECTIONS==cds.trials.tgtDir(trialNum)))) && ...
                    ~isempty(find(optsTask.BUMP_DIRECTIONS==cds.trials.bumpDir(trialNum))) && ~isempty(find(optsTask.STIM_CODE==cds.trials.stimCode(trialNum))))
                
                cellIdx = [find(strcmpi(optsTask.TRIAL_LIST,taskString)),find(optsTask.TARGET_DIRECTIONS==cds.trials.tgtDir(trialNum)),...
                    find(optsTask.BUMP_DIRECTIONS==cds.trials.bumpDir(trialNum)),find(optsTask.STIM_CODE==cds.trials.stimCode(trialNum))];
                
                % update counter for trial number
                trialCounter(cellIdx(1),cellIdx(2),cellIdx(3),cellIdx(4)) = trialCounter(cellIdx(1),cellIdx(2),cellIdx(3),cellIdx(4)) + 1;
                % get and store spike data
                spikeMask = cds.units(neuronNumber).spikes.ts > cds.trials.startTime(trialNum)-3 & cds.units(neuronNumber).spikes.ts < cds.trials.endTime(trialNum)+3;
                spikeTimes = cds.units(neuronNumber).spikes.ts(spikeMask) - cds.trials.(optsTask.ZERO_MARKER)(trialNum);
                if(~isempty(spikeTimes))
                    spikeTimeData{cellIdx(1),cellIdx(2),cellIdx(3),cellIdx(4)}(end+1:end+numel(spikeTimes),:) = ...
                        [spikeTimes];
                end  
                binMin = min(binMin,min(spikeTimes));
                binMax = max(binMax,max(spikeTimes));
                
                % populate stim times if available
                if(cds.trials.stimCode(trialNum) ~= -1) 
                    % for each stim on time, if its within trial start and end
                    % time, store it

                    for stimIdx = stimIdxStart:min(numel(cds.stimOn),stimIdxEnd)
                        if(cds.stimOn(stimIdx) > cds.trials.startTime(trialNum) && cds.stimOn(stimIdx) < cds.trials.endTime(trialNum))
                            stimData{cellIdx(1),cellIdx(2),cellIdx(3),cellIdx(4)}(end+1,1) = cds.stimOn(stimIdx) - cds.trials.(optsTask.ZERO_MARKER)(trialNum);
                            stimData{cellIdx(1),cellIdx(2),cellIdx(3),cellIdx(4)}(end,2) = trialCounter(cellIdx(1),cellIdx(2),cellIdx(3),cellIdx(4));
                            stimIdxStart = stimIdx;
                            stimIdxEnd = stimIdxStart + 100;
                        end
                    end
                    
                end
            end
        end
    end
    
    %% remove combinations that had too few trials
    for trialList = 1:numel(optsTask.TRIAL_LIST)
        for tgtDir = 1:numel(optsTask.TARGET_DIRECTIONS)
            for bumpDir = 1:numel(optsTask.BUMP_DIRECTIONS)
                for stimCode = 1:numel(optsTask.STIM_CODE)
                    if(trialCounter(trialList,tgtDir,bumpDir,stimCode) < 3)
                        trialCounter(trialList,tgtDir,bumpDir,stimCode) = 0;
                        spikeTimeData{trialList,tgtDir,bumpDir,stimCode} = [];
                        stimData{trialList,tgtDir,bumpDir,stimCode} = [];   
                    end
                end
            end
        end
    end
    %% Bin data
    binMin = binMin - mod(binMin,optsTask.BIN_SIZE) - 4*optsTask.BIN_SIZE;
    binMax = binMax - mod(binMax,optsTask.BIN_SIZE) + 4*optsTask.BIN_SIZE;
    % make sure that 0 is an edge of a bin
    bins = binMin:optsTask.BIN_SIZE:binMax;
    for trialList = 1:numel(optsTask.TRIAL_LIST)
        for tgtDir = 1:numel(optsTask.TARGET_DIRECTIONS)
            for bumpDir = 1:numel(optsTask.BUMP_DIRECTIONS)
                for stimCode = 1:numel(optsTask.STIM_CODE)
                    if(~isempty(spikeTimeData{trialList,tgtDir,bumpDir,stimCode})) % make sure to normalize by # trials
                        spikeTimeData{trialList,tgtDir,bumpDir,stimCode} = histcounts(spikeTimeData{trialList,tgtDir,bumpDir,stimCode},bins)/trialCounter(trialList,tgtDir,bumpDir,stimCode)/optsTask.BIN_SIZE;
                    end
                end
            end
        end
    end
    
    %% Combine data from same bumpDir or tgtDir according to optsTask.COMBINE
    for c = 1:numel(optsTask.COMBINE)
        switch optsTask.COMBINE{c}
            case 'tgtDir'
                spikeTimeDataTemp = cell(numel(optsTask.TRIAL_LIST),1,numel(optsTask.BUMP_DIRECTIONS),numel(optsTask.STIM_CODE));
                stimDataTemp = cell(numel(optsTask.TRIAL_LIST),1,numel(optsTask.BUMP_DIRECTIONS),numel(optsTask.STIM_CODE));
                
                for trialList = 1:numel(optsTask.TRIAL_LIST)
                    for bumpDir = 1:numel(optsTask.BUMP_DIRECTIONS)
                        for stimCode = 1:numel(optsTask.STIM_CODE)
                            for tgtDir = 1:numel(optsTask.TARGET_DIRECTIONS)
                                if(~isempty(spikeTimeData{trialList,tgtDir,bumpDir,stimCode}))
                                    if(~isempty(spikeTimeDataTemp{trialList,1,bumpDir,stimCode}))
                                        spikeTimeDataTemp{trialList,1,bumpDir,stimCode} = [spikeTimeDataTemp{trialList,1,bumpDir,stimCode};...
                                            spikeTimeData{trialList,tgtDir,bumpDir,stimCode}];
                                    else
                                         spikeTimeDataTemp{trialList,1,bumpDir,stimCode} = spikeTimeData{trialList,tgtDir,bumpDir,stimCode};
                                    end
                                    
                                    if(~isempty(stimDataTemp{trialList,1,bumpDir,stimCode}))
                                        stimDataTemp{trialList,1,bumpDir,stimCode} = [stimDataTemp{trialList,1,bumpDir,stimCode};...
                                            stimData{trialList,tgtDir,bumpDir,stimCode}];
                                    else
                                        stimDataTemp{trialList,1,bumpDir,stimCode} = stimData{trialList,tgtDir,bumpDir,stimCode};
                                    end
                                end
                            end
                        end
                    end
                end

                spikeTimeData = spikeTimeDataTemp;
                stimData = stimDataTemp;
                clear spikeTimeDataTemp;
                clear stimDataTemp;
                
            case 'bumpDir'
                spikeTimeDataTemp = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),1,numel(optsTask.STIM_CODE));
                stimDataTemp = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),1,numel(optsTask.STIM_CODE));
                
                for trialList = 1:numel(optsTask.TRIAL_LIST)
                    for tgtDir = 1:numel(optsTask.TARGET_DIRECTIONS)
                        for stimCode = 1:numel(optsTask.STIM_CODE)
                            for bumpDir = 1:numel(optsTask.BUMP_DIRECTIONS)
                                if(~isempty(spikeTimeData{trialList,tgtDir,bumpDir,stimCode}))
                                    if(~isempty(spikeTimeDataTemp{trialList,tgtDir,1,stimCode}))
                                        spikeTimeDataTemp{trialList,tgtDir,1,stimCode} = [spikeTimeDataTemp{trialList,tgtDir,1,stimCode};...
                                            spikeTimeData{trialList,tgtDir,bumpDir,stimCode}];
                                    else
                                         spikeTimeDataTemp{trialList,tgtDir,1,stimCode} = spikeTimeData{trialList,tgtDir,bumpDir,stimCode};
                                    end
                                    if(~isempty(stimDataTemp{trialList,tgtDir,1,stimCode}))
                                        stimDataTemp{trialList,tgtDir,1,stimCode} = [stimDataTemp{trialList,tgtDir,1,stimCode};...
                                            stimData{trialList,tgtDir,bumpDir,stimCode}];
                                    else
                                        stimDataTemp{trialList,tgtDir,1,stimCode} = stimData{trialList,tgtDir,bumpDir,stimCode};
                                    end
                                end
                            end
                        end
                    end
                end

                spikeTimeData = spikeTimeDataTemp;
                stimData = stimDataTemp;
                clear spikeTimeDataTemp
                clear stimDataTemp;
            case 'stimCode'
                spikeTimeDataTemp = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),numel(optsTask.BUMP_DIRECTIONS),1);
                stimDataTemp = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),numel(optsTask.BUMP_DIRECTIONS),1);
                
                for trialList = 1:numel(optsTask.TRIAL_LIST)
                    for tgtDir = 1:numel(optsTask.TARGET_DIRECTIONS)
                        for bumpDir = 1:numel(optsTask.BUMP_DIRECTIONS)
                            for stimCode = 1:numel(optsTask.STIM_CODE)
                                if(~isempty(spikeTimeData{trialList,tgtDir,bumpDir,stimCode}))
                                    if(~isempty(spikeTimeDataTemp{trialList,tgtDir,bumpDir,1}))
                                        spikeTimeDataTemp{trialList,tgtDir,bumpDir,1} = [spikeTimeDataTemp{trialList,tgtDir,bumpDir,1};...
                                            spikeTimeData{trialList,tgtDir,bumpDir,stimCode}];
                                    else
                                         spikeTimeDataTemp{trialList,tgtDir,bumpDir,1} = spikeTimeData{trialList,tgtDir,bumpDir,stimCode}(:,:);
                                    end
                                    if(~isempty(stimDataTemp{trialList,tgtDir,bumpDir,1}))
                                        stimDataTemp{trialList,tgtDir,bumpDir,1} = [stimDataTemp{trialList,tgtDir,bumpDir,1};...
                                            stimData{trialList,tgtDir,bumpDir,stimCode}];
                                    elseif(~isempty(stimData{trialList,tgtDir,bumpDir,stimCode}))
                                        stimDataTemp{trialList,tgtDir,bumpDir,1} = stimData{trialList,tgtDir,bumpDir,stimCode};
                                    end
                                end
                            end
                        end
                    end
                end
                spikeTimeData = spikeTimeDataTemp;
                stimData = stimDataTemp;
                clear spikeTimeDataTemp
                clear stimDataTemp;
        end
    end 
    
    %% now that we have the spike times, we can make rasters (yay)
    for trial = 1:size(spikeTimeData,1)
        for tgtDir = 1:size(spikeTimeData,2)
            for bumpDir = 1:size(spikeTimeData,3)
                for stimCode = 1:size(spikeTimeData,4)
                    % check to make sure there are trials
                    if(~isempty(spikeTimeData{trial,tgtDir,bumpDir,stimCode}))
                        if(optsPlot.PLOT_TITLE)
                            optsPlot.TITLE = strcat(optsTask.TRIAL_LIST{trial},' Target: ',num2str(optsTask.TARGET_DIRECTIONS(tgtDir)),...
                                ' Bump: ',num2str(optsTask.BUMP_DIRECTIONS(bumpDir)),' Stim: ',num2str(optsTask.STIM_CODE(stimCode)));
                        end
                        
                        if(strcmpi(optsTask.BIN_SIZE,'')==1)
                            optsTask.BIN_SIZE = 1/1000;
                        end
                        
                        %%%%%%%
                        %% fix this so that 0 is an edge of a bin
                        xData = bins(1:end-1) + (bins(2)-bins(1))/2;
                        yData = spikeTimeData{trial,tgtDir,bumpDir,stimCode};
                        xData = repmat(xData,size(yData,1),1);
                        
                        % if x limit is not specified, come up with one
                        % based on task start and end time
                        if(~optsPlot.xLimitFlag)
                            trialMask = cds.trials.result == 'R' & ~isnan(cds.trials.tgtOnTime) & ~isnan(cds.trials.movePeriod) & ~isnan(cds.trials.(optsTask.ZERO_MARKER));
                            trialLengths = cds.trials.endTime(trialMask) - cds.trials.startTime(trialMask);
                            [~,minTrialIdx] = min(trialLengths);
                            temp = find(trialMask); minTrialIdx = temp(minTrialIdx);
                            xlimits(1,1) = cds.trials.startTime(minTrialIdx) - cds.trials.(optsTask.ZERO_MARKER)(minTrialIdx);
                            xlimits(1,2) = cds.trials.endTime(minTrialIdx) - cds.trials.(optsTask.ZERO_MARKER)(minTrialIdx);
                            % round xlimits and expand a bit
                            xlimits = round(xlimits,1);
                            xlimits = [xlimits(1,1)-0.2,xlimits(1,2)+0.2];
                            optsPlot.X_LIMITS = xlimits;
                        end
                        % if y limit is not specified, come up with one
                        % based on the dimensions of yData
%                         ylimits(1,1) = 0;
%                         ylimits(1,2) = max(yData) + 0.8;
%                         optsPlot.Y_LIMITS = ylimits;
                        optsPlot.NUM_PLOTS = size(yData,1);
                        figureHandle{end+1} = plotPSTHLIB(xData',yData',optsPlot,optsSave);
                    end
                end
            end
        end
    end
    
    
    %% if optsTask.SAME_Y_LIMITS, then set all figure handles to the same y limits
    if(optsTask.SAME_Y_LIMITS)
        maxYLim = -10000;
        for fig = 1:numel(figureHandle)
            maxYLim = max(maxYLim,figureHandle{fig}.CurrentAxes.YLim(2));
        end
        for fig = 1:numel(figureHandle)
            figureHandle{fig}.CurrentAxes.YLim = [0,maxYLim];
        end
    end
    %% if optsTask.SAME_X_LIMITS, then set all figure handles to the same x limits
    if(optsTask.SAME_X_LIMITS)
        maxXLim = -10000;
        minXLim = 10000;
        for fig = 1:numel(figureHandle)
            maxXLim = max(maxXLim,figureHandle{fig}.CurrentAxes.XLim(2));
            minXLim = min(minXLim,figureHandle{fig}.CurrentAxes.XLim(1));
        end
        for fig = 1:numel(figureHandle)
            figureHandle{fig}.CurrentAxes.XLim = [minXLim,maxXLim];
        end
    end
end


function [optsPlot] = configureOptionsPlot(optsPlotInput)
    
    %% initialize possible fieldname variables
    optsPlot.MAKE_FIGURE = 1;
    optsPlot.X_LABEL = '';
    optsPlot.Y_LABEL = '';
    optsPlot.X_LIMITS = '';
    optsPlot.Y_LIMITS = '';
    optsPlot.X_TICK = '';
    optsPlot.Y_TICK = '';
    optsPlot.X_MINOR_TICK = '';
    optsPlot.Y_MINOR_TICK = '';
    optsPlot.X_TICK_LABEL = '';
    optsPlot.Y_TICK_LABEL = '';
    optsPlot.BAR_STYLE = 'bar';
    optsPlot.LINE_STYLE = '-';
    optsPlot.TITLE = '';
    optsPlot.MARKER_STYLE = '.';
    optsPlot.MARKER_COLOR = 'k';
    optsPlot.MARKER_SIZE = 4;
    optsPlot.LEGEND_STRING = '';
    optsPlot.PLOT_TITLE = 1;
    
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
    
    if(strcmpi(optsPlot.X_LIMITS,'')~=1)
        optsPlot.xLimitFlag = 1;
    else
        optsPlot.xLimitFlag = 0;
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


function [optsTask] = configureOptionsTask(optsTaskInput,cds)

    % default is to seperate by all bump directions and trial types
    optsTask.TRIAL_LIST = {}; 
    optsTask.TARGET_DIRECTIONS = unique(cds.trials.tgtDir(~isnan(cds.trials.tgtDir)));
    optsTask.BUMP_DIRECTIONS = [-1;unique(cds.trials.bumpDir(~isnan(cds.trials.bumpDir)))];
    optsTask.STIM_CODE = [-1;unique(cds.trials.stimCode(~isnan(cds.trials.stimCode)))];
    optsTask.ZERO_MARKER = 'goCueTime'; % needs to be a field in cds.trials
    optsTask.COMBINE = {'bumpDir'};
    optsTask.IGNORE = {};
    optsTask.PLOT_STIM_TIME = 0;
    optsTask.BIN_SIZE = '';
    optsTask.SAME_Y_LIMITS = 0;
    optsTask.SAME_X_LIMITS = 0;
    
    if(sum(cds.trials.ctrHoldBump==1 & isnan(cds.trials.stimCode)) > 0)
        optsTask.TRIAL_LIST{end+1,1} = 'ctrHoldBump';
    end
    if(sum(cds.trials.ctrHoldBump==1 & ~isnan(cds.trials.stimCode)) > 0)
        optsTask.TRIAL_LIST{end+1,1} = 'ctrHoldBumpStim';
    end
    if(sum(cds.trials.delayBump==1 & isnan(cds.trials.stimCode)) > 0)
        optsTask.TRIAL_LIST{end+1,1} = 'delayBump';
    end
    if(sum(cds.trials.delayBump==1 & ~isnan(cds.trials.stimCode)) > 0)
        optsTask.TRIAL_LIST{end+1,1} = 'delayBumpStim';
    end
    if(sum(cds.trials.moveBump==1 & isnan(cds.trials.stimCode)) > 0)
        optsTask.TRIAL_LIST{end+1,1} = 'moveBump';
    end
    if(sum(cds.trials.moveBump==1 & ~isnan(cds.trials.stimCode)) > 0)
        optsTask.TRIAL_LIST{end+1,1} = 'moveBumpStim';
    end
    
    %% check if in optsSave and optsSaveInput, overwrite if so
    try
        inputFieldnames = fieldnames(optsTaskInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(optsTask,inputFieldnames{fn}))
               optsTask.(inputFieldnames{fn}) = optsTaskInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
    
    %% parse through each field and remove things that don't exist?
end