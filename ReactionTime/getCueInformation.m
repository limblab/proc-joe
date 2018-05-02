function [cueInfo] = getCueInformation(cds,opts)

    % this function finds all of the cues in the data for a reaction time
    % task. Currently checks for:
    % stim sync line (assumes its ainp16)
    % motor control line/bump time word (these are at the same time, bump
    %   word is easier since its already been done)
    
    % this function then organizes that information into a useful struct
    % that is related to the cds trials table.
    % each entry into cue information will contain 
    % the trial idx in trials table
    % the cue onset time
    % any meta data about the cue (stim code, bump force, rise time, etc)
    
    %% configure opts and set default values
    opts = configureOpts(opts);
    cueInfo = [];
    
    %% find the stim sync line if it exists -- only supports ainp16 currently
    stimSyncIdx = [];
    stimSyncName = '';
    
    if(isempty(opts.STIM_SYNC_NAME) && isempty(opts.STIM_SYNC_ANALOG_IDX))
        %look for ainp16
        for j=1:numel(cds.analog)
            stimSyncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'ainp16'));
            if ~isempty(stimSyncIdx)
                stimSyncIdx=j;
                stimSyncName='ainp16';
            end
        end
    end
    
    %% get stim on time
    if(~isempty(stimSyncIdx))
        stimOn=cds.analog{stimSyncIdx}.t(find(diff(cds.analog{stimSyncIdx}.(stimSyncName)-mean(cds.analog{stimSyncIdx}.(stimSyncName))>3)>.5));
    end
    

    %% for each trial, enter into cueInfo
    cueInfo.cueTrueTime = NaN(size(cds.trials,1),1);
    cueInfo.cueTrialTime = NaN(size(cds.trials,1),1);
    cueInfo.stimCode = NaN(size(cds.trials,1),1);
    cueInfo.bumpMag = NaN(size(cds.trials,1),1);
    cueInfo.bumpRamp = NaN(size(cds.trials,1),1);
    cueInfo.bumpDir = NaN(size(cds.trials,1),1);
    
    for tr = 1:size(cds.trials,1)        
        
        % get cue time based on bump or stim
        % check for bump then check for stim
        if(~isnan(cds.trials.bumpTime(tr)))
            cueInfo.cueTrueTime(tr) = cds.trials.bumpTime(tr);
            cueInfo.cueTrialTime(tr) = cueInfo.cueTrueTime(tr) - cds.trials.(opts.START_VAR)(tr);
            cueInfo.bumpMag(tr) = cds.trials.bumpMagnitude(tr);
            cueInfo.bumpRamp(tr) = cds.trials.bumpRisePeriod(tr);
            cueInfo.bumpDir(tr) = cds.trials.bumpDir(tr);
        elseif(cds.trials.isStimTrial(tr) && ~isnan(cds.trials.stimCode(tr)) && ~isempty(find(stimOn > cds.trials.startTime(tr) & stimOn < cds.trials.endTime)))
            cueInfo.cueTrueTime(tr) = find(stimOn > cds.trials.startTime(tr) & stimOn < cds.trials.endTime,'first');
            cueInfo.cueTrialTime(tr) = cueInfo.cueTrueTime(tr) - cds.trials.(opts.START_VAR)(tr);
            cueInfo.stimCode(tr) = cds.trials.stimCode(tr);
        end        
        
                
    end
    

    
end



function [opts] = configureOpts(optsInput)

    opts = [];
    
    opts.STIM_SYNC_NAME = [];
    opts.STIM_SYNC_ANALOG_IDX = [];
    
    opts.START_VAR = 'tgtOnTime';
    %% check if in optsSave and optsSaveInput, overwrite if so
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