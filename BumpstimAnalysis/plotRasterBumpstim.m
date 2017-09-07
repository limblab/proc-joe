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
    %% extract data
    spikeTimeData = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),numel(optsTask.BUMP_DIRECTIONS),numel(optsTask.STIM_CODE));
    trialCounter = zeros(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),numel(optsTask.BUMP_DIRECTIONS),numel(optsTask.STIM_CODE));
    
    for trialNum = 1:size(cds.trials,1)
        if(cds.trials.result(trialNum) == 'R' && ~isnan(cds.trials.tgtOnTime(trialNum)) && ~isnan(cds.trials.movePeriod(trialNum))) % if a reward trial and not a corrupted trial
            taskString = '';
            
            if(cds.trials.ctrHoldBump(trialNum) && cds.trials.stimCode(trialNum)==-1 && sum(strcmpi(optsTask.TRIAL_LIST,'ctrHoldBump')) > 0)
                taskString = 'ctrHoldBump';
            elseif(cds.trials.ctrHoldBump(trialNum) && cds.trials.stimCode(trialNum)~=-1 && sum(strcmpi(optsTask.TRIAL_LIST,'ctrHoldBumpStim')) > 0)
                taskString = 'ctrHoldBumpStim';
            elseif(cds.trials.delayBump(trialNum) && cds.trials.stimCode(trialNum)==-1 && sum(strcmpi(optsTask.TRIAL_LIST,'delayBump')) > 0)
                taskString = 'delayBump';
            elseif(cds.trials.delayBump(trialNum) && cds.trials.stimCode(trialNum)~=-1 && sum(strcmpi(optsTask.TRIAL_LIST,'delayBumpStim')) > 0)
                taskString = 'delayBumpStim';
            elseif(cds.trials.moveBump(trialNum) && cds.trials.stimCode(trialNum)==-1 && sum(strcmpi(optsTask.TRIAL_LIST,'moveBump')) > 0)
                taskString = 'moveBump';
            elseif(cds.trials.moveBump(trialNum) && cds.trials.stimCode(trialNum)~=-1 && sum(strcmpi(optsTask.TRIAL_LIST,'moveBumpStim')) > 0)
                taskString = 'moveBumpStim';
            elseif(~cds.trials.ctrHoldBump(trialNum) && ~cds.trials.delayBump(trialNum) && ~cds.trials.moveBump(trialNum))
                taskString = 'noBump';
            end
            
            if(strcmpi(taskString,'')~=1 && ~isempty(find(optsTask.TARGET_DIRECTIONS==cds.trials.tgtDir(trialNum))) && ...
                    ~isempty(find(optsTask.BUMP_DIRECTIONS==cds.trials.bumpDir(trialNum))) && ~isempty(find(optsTask.STIM_CODE==cds.trials.stimCode(trialNum))))
                
                cellIdx = [find(strcmpi(optsTask.TRIAL_LIST,taskString)),find(optsTask.TARGET_DIRECTIONS==cds.trials.tgtDir(trialNum)),...
                    find(optsTask.BUMP_DIRECTIONS==cds.trials.bumpDir(trialNum)),find(optsTask.STIM_CODE==cds.trials.stimCode(trialNum))];
                
                % update counter for trial number
                trialCounter(cellIdx(1),cellIdx(2),cellIdx(3),cellIdx(4)) = trialCounter(cellIdx(1),cellIdx(2),cellIdx(3),cellIdx(4)) + 1;
                % get and store spike data
                spikeMask = cds.units(neuronNumber).spikes.ts > cds.trials.startTime(trialNum) & cds.units(neuronNumber).spikes.ts < cds.trials.endTime(trialNum);
                spikeTimes = cds.units(neuronNumber).spikes.ts(spikeMask) - cds.trials.(optsTask.ZERO_MARKER)(trialNum);
                if(~isempty(spikeTimes))
                    spikeTimeData{cellIdx(1),cellIdx(2),cellIdx(3),cellIdx(4)}(end+1:end+numel(spikeTimes),:) = ...
                        [spikeTimes,trialCounter(cellIdx(1),cellIdx(2),cellIdx(3),cellIdx(4))*ones(numel(spikeTimes),1)];
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
                    end
                end
            end
        end
     end
    
    %% Combine data from same bumpDir or tgtDir according to optsTask.COMBINE
    combineData = {};
    switch optsTask.COMBINE
        case 'tgtDir'
            spikeTimeDataTemp = cell(numel(optsTask.TRIAL_LIST),1,numel(optsTask.BUMP_DIRECTIONS),numel(optsTask.STIM_CODE));
            
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
                            end
                        end
                    end
                end
            end
            
            spikeTimeData = spikeTimeDataTemp;
            clear spikeTimeDataTemp;
        case 'bumpDir'
            spikeTimeDataTemp = cell(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),1,numel(optsTask.STIM_CODE));
            trialCounterTemp = zeros(numel(optsTask.TRIAL_LIST),numel(optsTask.TARGET_DIRECTIONS),1,numel(optsTask.STIM_CODE));
            
            for trialList = 1:numel(optsTask.TRIAL_LIST)
                for tgtDir = 1:numel(optsTask.TARGET_DIRECTIONS)
                    for stimCode = 1:numel(optsTask.STIM_CODE)
                        for bumpDir = 1:numel(optsTask.BUMP_DIRECTIONS)
                            if(~isempty(spikeTimeData{trialList,tgtDir,bumpDir,stimCode}))
                                if(~isempty(spikeTimeDataTemp{trialList,tgtDir,1,stimCode}))
                                    combineData{trialList,tgtDir,1,stimCode}(end+1,1) = spikeTimeDataTemp{trialList,1,bumpDir,stimCode}(end,2);
                                    spikeTimeDataTemp{trialList,tgtDir,1,stimCode} = [spikeTimeDataTemp{trialList,tgtDir,1,stimCode};...
                                        spikeTimeData{trialList,tgtDir,bumpDir,stimCode}(:,1),spikeTimeData{trialList,tgtDir,bumpDir,stimCode}(:,2)+spikeTimeDataTemp{trialList,tgtDir,1,stimCode}(end,2)];
                                else
                                     spikeTimeDataTemp{trialList,1,bumpDir,stimCode} = spikeTimeData{trialList,tgtDir,bumpDir,stimCode}(:,:);
                                end
                            end
                        end
                    end
                end
            end
            
            spikeTimeData = spikeTimeDataTemp;
            clear spikeTimeDataTemp
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
                            optsPlot.DIVIDING_LINES = combineData{trial,tgtDir,bumpDir,stimCode};
                            optsPlot.DIVIDING_LINES_COLORS = {'r','b','g','m','k'};
                        end
                        xData = spikeTimeData{trial,tgtDir,bumpDir,stimCode}(:,1);
                        yData = spikeTimeData{trial,tgtDir,bumpDir,stimCode}(:,2);
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
    optsTask.COMBINE = 'tgtDir';
    
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
        inputFieldnames = fieldnames(optsBumpInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(optsTask,inputFieldnames{fn}))
               optsTask.(inputFieldnames{fn}) = optsBumpInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
    
    %% parse through each field and remove things that don't exist?
end