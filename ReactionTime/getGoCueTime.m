function [td] = getGoCueTime(td,cds)

    % this function finds the go cue time and puts an entry into trial data
    % that is the true go cue time and the corresponding idx. This is used
    % for reaction time as the go cues from the behavior data do not align
    % with the actual cue presentation. Instead, the analog sync lines are
    % used for the stim input. The bump time is used for the bump cue time.

    
    %% find the stim and visual sync line if it exists -- only supports ainp16 currently
    stimSyncIdx = [];
    stimSyncName = '';
    
    visualSyncIdx = [];
    visualSyncName = '';
    %look for ainp16 for stim sync or videosync....
    for j=1:numel(cds.analog)
        stimSyncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'ainp16'));
        if ~isempty(stimSyncIdx)
            stimSyncIdx=j;
            stimSyncName='ainp16';
        end
    end
    
    if(isempty(stimSyncIdx))
        for j=1:numel(cds.analog)
            stimSyncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'videosync'));
            if ~isempty(stimSyncIdx)
                stimSyncIdx=j;
                stimSyncName='videosync';
            end
        end
    end
    
    %look for visualTgtSync for visual sync
    for j=1:numel(cds.analog)
        visualSyncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'visualTgtSync'));
        if ~isempty(visualSyncIdx)
            visualSyncIdx=j;
            visualSyncName='visualTgtSync';
        end
    end
    
    if(isempty(visualSyncIdx))
        for j=1:numel(cds.analog)
            visualSyncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'viusalTgt'));
            if ~isempty(visualSyncIdx)
                visualSyncIdx=j;
                visualSyncName='viusalTgt';
            end
        end
    end
    
    if(isempty(visualSyncIdx))
        for j=1:numel(cds.analog)
            visualSyncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'visualsync'));
            if ~isempty(visualSyncIdx)
                visualSyncIdx=j;
                visualSyncName='visualsync';
            end
        end
    end
    
    if(isempty(visualSyncIdx))
        for j=1:numel(cds.analog)
            visualSyncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'ainp15'));
            if ~isempty(visualSyncIdx)
                visualSyncIdx=j;
                visualSyncName='ainp15';
            end
        end
    end
    
    %% get stim on time
    stimOn = [];
    if(~isempty(stimSyncIdx))
        stimOn=cds.analog{stimSyncIdx}.t(find(diff(cds.analog{stimSyncIdx}.(stimSyncName)-mean(cds.analog{stimSyncIdx}.(stimSyncName))>3)>.5));
    end
    
    %% get visual on time
%     visualOn = [];
%     if(~isempty(visualSyncIdx))
%         visualOn=cds.analog{visualSyncIdx}.t(find(diff(cds.analog{visualSyncIdx}.(visualSyncName)-mean(cds.analog{visualSyncIdx}.(visualSyncName))>10)>0.5));
%     end

    %% set idx_goCueTime and goCueTime in td
    for tr = 1:numel(td)
        % determine bump time for this trial
        bumpTime = cds.trials.bumpTime(td(tr).trial_id)-cds.trials.startTime(td(tr).trial_id);
        
        % determine stim on time for this trial
        stimTime = [];
        if(~isempty(stimSyncIdx))
            stimTime = stimOn(find(stimOn > cds.trials.startTime(td(tr).trial_id) & stimOn < cds.trials.endTime(td(tr).trial_id),1,'first'));
            stimTime = stimTime - cds.trials.startTime(td(tr).trial_id);
        end
        
        if(~isnan(bumpTime)) %prioritize bump time as go cue as stim has some delays
            td(tr).goCueTime = bumpTime + td(tr).idx_startTime*td(tr).bin_size;
            td(tr).idx_goCueTime = floor(bumpTime/td(tr).bin_size + td(tr).idx_startTime);
        elseif(td(tr).isStimTrial && ~isempty(stimTime))%check for stim only trial
            td(tr).goCueTime = stimTime + td(tr).idx_startTime*td(tr).bin_size;
            td(tr).idx_goCueTime = floor(stimTime/td(tr).bin_size + td(tr).idx_startTime);
            
        elseif(~isnan(td(tr).isVisualTrial) && td(tr).isVisualTrial)
            trialGoCueTime = cds.trials.goCueTime(td(tr).trial_id);
            % get visual onset time
            if(~isempty(visualSyncIdx))
                mean_sync = mean(cds.analog{visualSyncIdx}.(visualSyncName)([find(cds.analog{visualSyncIdx}.t > trialGoCueTime,1,'first'),...
                    find(cds.analog{visualSyncIdx}.t > trialGoCueTime + 0.04,1,'first')]));
                visualOn=cds.analog{visualSyncIdx}.t(find(diff(cds.analog{visualSyncIdx}.(visualSyncName)-mean_sync>10)>0.5));
            else
                visualOn = trialGoCueTime;
            end
            visualTime = visualOn(find(visualOn > trialGoCueTime,1,'first'));
            visualTime = visualTime - cds.trials.startTime(td(tr).trial_id) + 5/1000; % 8 ms is the measured time between when the corner tgt (the sync line) is on and when t
            %he center tgt actually turns on, 3 ms is the time after the tgt turns on that we detect it turns on (8-3=5)
            
            if(isempty(visualTime) || visualTime - trialGoCueTime > 1)
                td(tr).goCueTime = nan;
                td(tr).idx_goCueTime = nan;
            else
                td(tr).goCueTime = visualTime + td(tr).idx_startTime*td(tr).bin_size;
                td(tr).idx_goCueTime = floor(visualTime/td(tr).bin_size + td(tr).idx_startTime);
            end
        else % if bumpTime is NaN && stimTime is empty, then no cue was
             % presented, set idx_goCueTime and goCueTime to NaN   
            td(tr).idx_goCueTime = NaN;
            td(tr).goCueTime = NaN;
        end
        
        if(~isnan(td(tr).idx_goCueTime) && td(tr).idx_goCueTime > td(tr).idx_endTime) % sanitize data
            td(tr).idx_goCueTime = NaN;
            td(tr).goCueTime = NaN;
        end
    end


end