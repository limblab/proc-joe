function [ figureHandle ] = plotRasterBumpstim( cds,neuronNumber,optsTaskInput,optsPlotInput,optsSaveInput )
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
                        [spikeTimes,trialCounter(cellIdx(1),cellIdx(2),cellIdx(3),cellIdx(4))*ones(numel(spikeTimes),1)];
                end  
                
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
    
    %% Combine data from same bumpDir or tgtDir according to optsTask.COMBINE
    for c = 1:numel(optsTask.COMBINE)
        switch optsTask.COMBINE{c}
            case 'tgtDir'
                combineData = cell(numel(optsTask.TRIAL_LIST),1,numel(optsTask.BUMP_DIRECTIONS),numel(optsTask.STIM_CODE));
                spikeTimeDataTemp = cell(numel(optsTask.TRIAL_LIST),1,numel(optsTask.BUMP_DIRECTIONS),numel(optsTask.STIM_CODE));
                stimDataTemp = cell(numel(optsTask.TRIAL_LIST),1,numel(optsTask.BUMP_DIRECTIONS),numel(optsTask.STIM_CODE));
                
                for trialList = 1:numel(optsTask.TRIAL_LIST)
                    for bumpDir = 1:numel(optsTask.BUMP_DIRECTIONS)
                        for stimCode = 1:numel(optsTask.STIM_CODE)
                            for tgtDir = 1:numel(optsTask.TARGET_DIRECTIONS)
                                if(~isempty(spikeTimeData{trialList,tgtDir,bumpDir,stimCode}))
                                    if(~isempty(spikeTimeDataTemp{trialList,1,bumpDir,stimCode}))
                                        combineData{trialList,1,bumpDir,stimCode}(end+1,1) = spikeTimeDataTemp{trialList,1,bumpDir,stimCode}(end,2);
                                        spikeTimeDataTemp{trialList,1,bumpDir,stimCode} = [spikeTimeDataTemp{trialList,1,bumpDir,stimCode};...
                                            spikeTimeData{trialList,tgtDir,bumpDir,stimCode}(:,1),spikeTimeData{trialList,tgtDir,bumpDir,stimCode}(:,2)+spikeTimeDataTemp{trialList,1,bumpDir,stimCode}(end,2)];
                                    else
                                         spikeTimeDataTemp{trialList,1,bumpDir,stimCode} = spikeTimeData{trialList,tgtDir,bumpDir,stimCode}(:,:);
                                    end
                                    
                                    if(~isempty(stimDataTemp{trialList,1,bumpDir,stimCode}))
                                        stimDataTemp{trialList,1,bumpDir,stimCode} = [stimDataTemp{trialList,1,bumpDir,stimCode};...
                                            stimData{trialList,tgtDir,bumpDir,stimCode}(:,1),stimData{trialList,tgtDir,bumpDir,stimCode}(:,2)+stimDataTemp{trialList,1,bumpDir,stimCode}(end,2)];
                                    else
                                        stimDataTemp{trialList,1,bumpDir,stimCode} = stimData{trialList,tgtDir,bumpDir,stimCode}(:,:);
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
                combineData = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),1,numel(optsTask.STIM_CODE));
                spikeTimeDataTemp = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),1,numel(optsTask.STIM_CODE));
                stimDataTemp = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),1,numel(optsTask.STIM_CODE));
                
                for trialList = 1:numel(optsTask.TRIAL_LIST)
                    for tgtDir = 1:numel(optsTask.TARGET_DIRECTIONS)
                        for stimCode = 1:numel(optsTask.STIM_CODE)
                            for bumpDir = 1:numel(optsTask.BUMP_DIRECTIONS)
                                if(~isempty(spikeTimeData{trialList,tgtDir,bumpDir,stimCode}))
                                    if(~isempty(spikeTimeDataTemp{trialList,tgtDir,1,stimCode}))
                                        combineData{trialList,tgtDir,1,stimCode}(end+1,1) = spikeTimeDataTemp{trialList,tgtDir,1,stimCode}(end,2);
                                        spikeTimeDataTemp{trialList,tgtDir,1,stimCode} = [spikeTimeDataTemp{trialList,tgtDir,1,stimCode};...
                                            spikeTimeData{trialList,tgtDir,bumpDir,stimCode}(:,1),spikeTimeData{trialList,tgtDir,bumpDir,stimCode}(:,2)+spikeTimeDataTemp{trialList,tgtDir,1,stimCode}(end,2)];
                                    else
                                         spikeTimeDataTemp{trialList,tgtDir,1,stimCode} = spikeTimeData{trialList,tgtDir,bumpDir,stimCode}(:,:);
                                    end
                                    if(~isempty(stimDataTemp{trialList,tgtDir,1,stimCode}))
                                        stimDataTemp{trialList,tgtDir,1,stimCode} = [stimDataTemp{trialList,tgtDir,1,stimCode};...
                                            stimData{trialList,tgtDir,bumpDir,stimCode}(:,1),stimData{trialList,tgtDir,bumpDir,stimCode}(:,2)+stimDataTemp{trialList,tgtDir,1,stimCode}(end,2)];
                                    else
                                        stimDataTemp{trialList,tgtDir,1,stimCode} = stimData{trialList,tgtDir,bumpDir,stimCode}(:,:);
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
                combineData = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),numel(optsTask.BUMP_DIRECTIONS),1);
                combineData(:,:,:,:) = {0};
                spikeTimeDataTemp = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),numel(optsTask.BUMP_DIRECTIONS),1);
                stimDataTemp = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),numel(optsTask.BUMP_DIRECTIONS),1);
                
                for trialList = 1:numel(optsTask.TRIAL_LIST)
                    for tgtDir = 1:numel(optsTask.TARGET_DIRECTIONS)
                        for bumpDir = 1:numel(optsTask.BUMP_DIRECTIONS)
                            for stimCode = 1:numel(optsTask.STIM_CODE)
                                if(~isempty(spikeTimeData{trialList,tgtDir,bumpDir,stimCode}))
                                    if(~isempty(spikeTimeDataTemp{trialList,tgtDir,bumpDir,1}))
                                        combineData{trialList,tgtDir,bumpDir,1}(end+1,1) = spikeTimeDataTemp{trialList,tgtDir,bumpDir,1}(end,2);
                                        spikeTimeDataTemp{trialList,tgtDir,bumpDir,1} = [spikeTimeDataTemp{trialList,tgtDir,bumpDir,1};...
                                            spikeTimeData{trialList,tgtDir,bumpDir,stimCode}(:,1),spikeTimeData{trialList,tgtDir,bumpDir,stimCode}(:,2)+spikeTimeDataTemp{trialList,tgtDir,bumpDir,1}(end,2)];
                                    else
                                         spikeTimeDataTemp{trialList,tgtDir,bumpDir,1} = spikeTimeData{trialList,tgtDir,bumpDir,stimCode}(:,:);
                                    end
                                    if(~isempty(stimDataTemp{trialList,tgtDir,bumpDir,1}))
                                        stimDataTemp{trialList,tgtDir,bumpDir,1} = [stimDataTemp{trialList,tgtDir,bumpDir,1};...
                                            stimData{trialList,tgtDir,bumpDir,stimCode}(:,1),stimData{trialList,tgtDir,bumpDir,stimCode}(:,2)+stimDataTemp{trialList,tgtDir,bumpDir,1}(end,2)];
                                    elseif(~isempty(stimData{trialList,tgtDir,bumpDir,stimCode}))
                                        stimDataTemp{trialList,tgtDir,bumpDir,1} = [stimData{trialList,tgtDir,bumpDir,stimCode}(:,1),...
                                            stimData{trialList,tgtDir,bumpDir,stimCode}(:,2) + combineData{trialList,tgtDir,bumpDir,1}(end)];
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
%     spikeTimeData = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),numel(optsTask.BUMP_DIRECTIONS),numel(optsTask.STIM_CODE));
 
    %% now that we have the spike times, we can make rasters (yay)
    for trial = 1:size(spikeTimeData,1)
        for tgtDir = 1:size(spikeTimeData,2)
            for bumpDir = 1:size(spikeTimeData,3)
                for stimCode = 1:size(spikeTimeData,4)
                    % check to make sure there are trials
                    if(~isempty(spikeTimeData{trial,tgtDir,bumpDir,stimCode}))
                        %% fix the title text later -- joe
                        optsPlot.TITLE = strcat(optsTask.TRIAL_LIST{trial},' Target: ',num2str(optsTask.TARGET_DIRECTIONS(tgtDir)),...
                            ' Bump: ',num2str(optsTask.BUMP_DIRECTIONS(bumpDir)),' Stim: ',num2str(optsTask.STIM_CODE(stimCode)));
                        if(~isempty(optsTask.COMBINE))
                            % remove 0 entry if there
                            if(~isempty(combineData{trial,tgtDir,bumpDir,stimCode}) && combineData{trial,tgtDir,bumpDir,stimCode}(1) == 0)
                                combineData{trial,tgtDir,bumpDir,stimCode} = combineData{trial,tgtDir,bumpDir,stimCode}(2:end);
                            end
                            
                            optsPlot.DIVIDING_LINES = combineData{trial,tgtDir,bumpDir,stimCode};
                            optsPlot.DIVIDING_LINES_COLORS = {'r','b',[0 0.6 0],'m','k'};
                        end
                        xData = spikeTimeData{trial,tgtDir,bumpDir,stimCode}(:,1);
                        yData = spikeTimeData{trial,tgtDir,bumpDir,stimCode}(:,2);
                        
                        % if plot stim times, put
                        % stim times into optsPlot. Sort in time for each
                        % grouping
                        if(optsTask.PLOT_STIM_TIME && ~isempty(stimData{trial,tgtDir,bumpDir,stimCode}))
                            optsPlot.PLOT_STIM_TIME = 1;
                            optsPlot.STIM_DATA_X = stimData{trial,tgtDir,bumpDir,stimCode}(:,1);
                            optsPlot.STIM_DATA_Y = stimData{trial,tgtDir,bumpDir,stimCode}(:,2);
                            
                            % only plot 1st stimulation time
                            
                            [optsPlot.STIM_DATA_Y,xIdx,~] = unique(optsPlot.STIM_DATA_Y,'first');
                            optsPlot.STIM_DATA_X = optsPlot.STIM_DATA_X(xIdx);
                            
                            %% sort trials within a combined group by
                            if(~isempty(optsTask.COMBINE))
                                combineMarkers = combineData{trial,tgtDir,bumpDir,stimCode};
                            else
                                combineMarkers = [];
                            end
                            
                            % reformat stim data if necessary
                            if(~isempty(optsPlot.STIM_DATA_Y) && max(optsPlot.STIM_DATA_Y) ~= size(optsPlot.STIM_DATA_Y,1))
                                newStimDataY = [1:1:max(yData)]';
                                newStimDataX = [-10000*ones(max(yData),1)];
                                newStimDataX(optsPlot.STIM_DATA_Y) = optsPlot.STIM_DATA_X;
                                optsPlot.STIM_DATA_X = newStimDataX;
                                optsPlot.STIM_DATA_Y = newStimDataY;
                            end

                            % sort trials by stim time within a combined
                            % trial
                            for cMark = 1:numel(combineMarkers)+1
                                % initial stim time (stim_data_y)
                                if(cMark == 1 && isempty(combineMarkers))
                                    [optsPlot.STIM_DATA_X,sortIdx] = sort(optsPlot.STIM_DATA_X);
                                    offset = 0;
                                elseif(cMark == 1)
                                    [optsPlot.STIM_DATA_X(1:combineMarkers(1)-1),sortIdx] = sort(optsPlot.STIM_DATA_X(1:combineMarkers(1)-1));
                                    offset = 0;
                                elseif(cMark == numel(combineMarkers) + 1)
                                    [optsPlot.STIM_DATA_X(combineMarkers(cMark-1):end),sortIdx] = sort(optsPlot.STIM_DATA_X(combineMarkers(cMark-1):end));
                                    offset = combineMarkers(cMark-1);
                                else
                                    [optsPlot.STIM_DATA_X(combineMarkers(cMark-1):combineMarkers(cMark)-1),sortIdx] = sort(optsPlot.STIM_DATA_X(combineMarkers(cMark-1):combineMarkers(cMark)-1));
                                    offset = combineMarkers(cMark-1);
                                end
                                % sort yData according to sortIdx
                                yDataTemp = yData;
                                for idx = 1:max(sortIdx)
                                    yDataTemp(yData==idx+offset) = sortIdx(idx)+offset;
                                end
                                yData = yDataTemp;
                            end
                                    
                        end
                        
                        
                        % if x limit is not specified, come up with one
                        % based on task start and end time
                        if(~optsPlot.xLimitFlag)
                            trialMask = cds.trials.result == 'R' & ~isnan(cds.trials.tgtOnTime) & ~isnan(cds.trials.movePeriod) & ~isnan(cds.trials.(optsTask.ZERO_MARKER));
                            trialLengths = cds.trials.endTime(trialMask) - cds.trials.startTime(trialMask);
                            [~,minTrialIdx] = min(trialLengths);
                            temp = find(trialMask); minTrialIdx = temp(minTrialIdx);
                            xLimits(1,1) = cds.trials.startTime(minTrialIdx) - cds.trials.(optsTask.ZERO_MARKER)(minTrialIdx);
                            xLimits(1,2) = cds.trials.endTime(minTrialIdx) - cds.trials.(optsTask.ZERO_MARKER)(minTrialIdx);
                            % round xlimits and expand a bit
                            xLimits = round(xLimits,1);
                            xLimits = [xLimits(1,1)-0.3,xLimits(1,2)+0.3];
                            optsPlot.X_LIMITS = xLimits;
                        end
                        % remove data that is before or after the xlimits
                        yData = yData(xData > optsPlot.X_LIMITS(1) & xData < optsPlot.X_LIMITS(2));
                        xData = xData(xData > optsPlot.X_LIMITS(1) & xData < optsPlot.X_LIMITS(2));
                        
                        % if y limit is not specified, come up with one
                        % based on the dimensions of yData
                        if(~optsPlot.yLimitFlag)
                            yLimits(1,1) = -0.8;
                            yLimits(1,2) = max(yData) + 0.8;
                            optsPlot.Y_LIMITS = yLimits;
                        end
                        
                        if(~isempty(combineData{trial,tgtDir,bumpDir,stimCode}))
                            optsPlot.Y_TICK = [1,combineData{trial,tgtDir,bumpDir,stimCode}',max(yData)];
                        else
                            optsPlot.Y_TICK = [1,max(yData)];
                        end
                        optsPlot.Y_MINOR_TICK = 'off';
                        plotRasterLIB(xData,yData,optsPlot,optsSave);
                        
                    end
                end
            end
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
    optsPlot.LINE_STYLE = '';
    optsPlot.TITLE = '';
    optsPlot.MARKER_STYLE = '.';
    optsPlot.MARKER_COLOR = 'k';
    optsPlot.MARKER_SIZE = 4;
    
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
    if(strcmpi(optsPlot.Y_LIMITS,'')~=1)
        optsPlot.yLimitFlag = 1;
    else
        optsPlot.yLimitFlag = 0;
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
    if(sum(cds.trials.ctrHoldBump~=1 & cds.trials.delayBump ~= 1 & cds.trials.moveBump~=1) > 0)
        optsTask.TRIAL_LIST{end+1,1} = 'noBump';
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
    
    % format COMBINE
    if(~iscell(optsTask.COMBINE))
        optsTask.COMBINE = {optsTask.COMBINE};
    end
    
    % format IGNORE and remake other things if necessary
    if(~iscell(optsTask.IGNORE))
        optsTask.IGNORE = {optsTask.IGNORE};
    end
    
end