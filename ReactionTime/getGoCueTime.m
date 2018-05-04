function [td] = getGoCueTime(td,cds)

    % this function finds the go cue time and puts an entry into trial data
    % that is the true go cue time and the corresponding idx. This is used
    % for reaction time as the go cues from the behavior data do not align
    % with the actual cue presentation. Instead, the analog sync lines are
    % used for the stim input. The bump time is used for the bump cue time.

    
    %% find the stim sync line if it exists -- only supports ainp16 currently
    stimSyncIdx = [];
    stimSyncName = '';
    
    %look for ainp16
    for j=1:numel(cds.analog)
        stimSyncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'ainp16'));
        if ~isempty(stimSyncIdx)
            stimSyncIdx=j;
            stimSyncName='ainp16';
        end
    end
    
    %% get stim on time
    if(~isempty(stimSyncIdx))
        stimOn=cds.analog{stimSyncIdx}.t(find(diff(cds.analog{stimSyncIdx}.(stimSyncName)-mean(cds.analog{stimSyncIdx}.(stimSyncName))>3)>.5));
    end

    %% set idx_goCueTime and goCueTime in td
    for tr = 1:numel(td)
        % determine bump time for this trial
        bumpTime = cds.trials.bumpTime(td(tr).trial_id)-cds.trials.startTime(td(tr).trial_id);
        
        % determine stim on time for this trial
        stimTime = stimOn(find(stimOn > cds.trials.startTime(td(tr).trial_id) & stimOn < cds.trials.endTime(td(tr).trial_id),1,'first'));
        stimTime = stimTime - cds.trials.startTime(td(tr).trial_id);
        
        if(~isnan(bumpTime)) %prioritize bump time as go cue as stim has some delays
            td(tr).goCueTime = bumpTime + td(tr).idx_startTime*td(tr).bin_size;
            td(tr).idx_goCueTime = floor(bumpTime/td(tr).bin_size + td(tr).idx_startTime);
        elseif(td(tr).isStimTrial && ~isempty(stimTime))%check for stim only trial
            td(tr).goCueTime = stimTime + td(tr).idx_startTime*td(tr).bin_size;
            td(tr).idx_goCueTime = floor(stimTime/td(tr).bin_size + td(tr).idx_startTime);
        else % if bumpTime is NaN && stimTime is empty, then no cue was
             % presented, set idx_goCueTime and goCueTime to NaN   
            td(tr).idx_goCueTime = NaN;
            td(tr).goCueTime = NaN;
        end
    end


end