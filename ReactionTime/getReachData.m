function [reachData] = getReachData(cds,cueInfo,opts)

    % this function gets all of the reaction times, movement onset and
    % kinematics for rewarded trials with cues (info provided by cueInfo)

    %% configure opts
    opts = configureOpts(opts);
    reachData = [];
    reachData.preTime = opts.PRE_TIME;
    reachData.postTime = opts.POST_TIME;
    
    %% get reaches for each trial
    for tr = 1:size(cds.trials,1)
        % initialize reach data parameters
        reachData.kin(tr,1).x = [];
        reachData.kin(tr,1).y = [];
        reachData.kin(tr,1).vx = [];
        reachData.kin(tr,1).vy = [];
        reachData.kin(tr,1).ax = [];
        reachData.kin(tr,1).ay = [];
        reachData.kin(tr,1).moveOnIdx = [];
        reachData.kin(tr,1).t = [];
        reachData.reactionTime(tr,1) = NaN;
        reachData.goCueTime(tr,1) = NaN;
        
        % populate reach data parameters
        if(~isnan(cueInfo.cueTrialTime(tr)) && cds.trials.result(tr) == 'R') % if we have a cue time and it was a rewarded trial. Discarding trials where he doesn't react
            % find kin start and end idx
            startIdx = find(cds.kin.t > cds.trials.(opts.START_VAR)(tr)+opts.PRE_TIME,1,'first');
            endIdx = find(cds.kin.t > cds.trials.(opts.END_VAR)(tr)+opts.POST_TIME,1,'first');
            
            % populate kin data
            reachData.kin(tr).x = cds.kin.x(startIdx:endIdx);
            reachData.kin(tr).y = cds.kin.y(startIdx:endIdx);
            reachData.kin(tr).vx = cds.kin.vx(startIdx:endIdx);
            reachData.kin(tr).vy = cds.kin.vy(startIdx:endIdx);
            reachData.kin(tr).ax = cds.kin.ax(startIdx:endIdx);
            reachData.kin(tr).ay = cds.kin.ay(startIdx:endIdx);
            reachData.kin(tr).t = cds.kin.t(startIdx:endIdx) - cds.trials.(opts.START_VAR)(tr);
            
            % put cue time in relative to start_var
            reachData.goCueTime(tr) = cueInfo.cueTrueTime(tr) - cds.trials.(opts.START_VAR)(tr);
            
            % find movement on idx and time
            opts.TARGET_DIRECTION = cds.trials.tgtDir(tr);
            reachData.kin(tr).moveOnIdx = getMovementOnset(reachData.kin(tr),cueInfo.cueTrialTime(tr),opts);
            if(isnan(reachData.kin(tr).moveOnIdx))
                reachData.kin(tr).moveOnTime = NaN;
                disp('could not find a move on time')
            else
                reachData.kin(tr).moveOnTime = reachData.kin(tr).t(reachData.kin(tr).moveOnIdx);
                if(reachData.kin(tr).moveOnTime-cueInfo.cueTrialTime(tr) > opts.MAX_REACTION_TIME)
                    reachData.kin(tr).moveOnIdx = NaN;
                    reachData.kin(tr).moveOnTime = NaN;
                end
            end
            % find reaction time
            reachData.reactionTime(tr) = reachData.kin(tr).moveOnTime - reachData.goCueTime(tr);
            
        end
    end
    
    
    
    %% now break reach data up based on cue 
    % assumption -- bumpMag and stimCode will cover all of the cues
    
    if(isempty(opts.BUMP_MAGS))
        bumpMags = unique(cueInfo.bumpMag(~isnan(cueInfo.bumpMag)));
    else
        bumpMags = opts.BUMP_MAGS;
    end
    
    if(isempty(opts.STIM_CODES))
        stimCodes = unique(cueInfo.stimCode(~isnan(cueInfo.stimCode)));
    else
        stimCodes = opts.STIM_CODES;
    end
    
    % split reaction times up based on cue, store reaction time
    cueIdx = 1;
    reactionTimse = {};
    %  each bump magnitude
    for bM = 1:numel(bumpMags)
        reactionTimes{cueIdx} = reachData.reactionTime(find(isEqual(cueInfo.bumpMag,bumpMags(bM)) & ~isnan(reachData.reactionTime)));
        cueIdx = cueIdx + 1;
    end
    
    %  each stim code
    for sC = 1:numel(stimCodes)
        reactionTimes{cueIdx} = reachData.reactionTime(isEqual(cueInfo.stimCode,stimCodes(sC)&~isnan(reachData.reactionTime)));
        cueIdx = cueIdx + 1;
    end
    
    % store relavant data into reachData
    reachData.bumpMags = [bumpMags;NaN(numel(stimCodes),1)];
    reachData.stimCodes = [NaN(numel(bumpMags),1);stimCodes];
    reachData.reactionTimes = reactionTimes;
    
end

function [out] = isEqual(data1,data2,thresh)

    out = data1 < data2 + eps & data1 > data2 - eps;

end

function [moveOnIdx] = getMovementOnset(kin,cueTime,opts)
    
    % takes in kinematic information and options
    % outputs when movement is on.
    
    % options:
    %   .MOVE_START_OFFSET = indices after cueTime at which a movement onset
    %   could be
    %   .TARGET_DIRECTION = target direction (in deg). If set to anything useful,
    %   then we only look at movements in that direction --- CURRENTLY NOT
    %   IMPLEMENTED, ASSUMES TARGET is TO THE RIGHT (0 deg)
    moveOnIdx = NaN; % not a number if a move on time could not be found

    % moveOnMask representes potential movement onset positions
    moveOnMask = ones(numel(kin.x),1);
    % remove those before the cue and offset
    cueIdx = find(kin.t > cueTime,1,'first');
    moveOnMask(1:(cueIdx+opts.MOVE_START_OFFSET)) = 0;
    % with remaining data, find peak vel baesd on accel crossing zero,
    % finds first peak
    s = kin.vx;
    ds = [0;diff(s)];
    dds = [0;diff(ds)];
    switch opts.METHOD
        case 'peakAcceleration'
%             peaks = [dds(1:end-1) > 0 & dds(2:end) < 0; 0]; % find peak acceleration
%             peakIdx = find(peaks & ds > opts.MIN_VAR & moveOnMask,1,'first'); % peak acceleration is first peak in the data range
            [~,peakIdx] = max(ds(cueIdx:cueIdx+opts.MOVE_END_OFFSET)); % find peak acceleration, make sure its within a certain range
            
            if(~isempty(peakIdx))
                peakIdx = peakIdx + cueIdx;
                % find latest time at which the acceleration is half maximum --
                % define that as movement onset
                thresh = ds(peakIdx)*opts.THRESH_MULT;
                moveOnIdx = find(moveOnMask & ds < thresh & (1:numel(moveOnMask))' < peakIdx,1,'last');

            end
        case 'peakVelocity'
            peaks = [ds(1:end-1) > 0 & ds(2:end) < 0; 0]; % find peak acceleration
            peakIdx = find(peaks & s > opts.MIN_VAR & moveOnMask,1,'first'); % peak acceleration is first peak in the data range

            if(~isempty(peakIdx))
                % find latest time at which the acceleration is half maximum --
                % define that as movement onset
                thresh = s(peakIdx)*opts.THRESH_MULT;
                moveOnIdx = find(moveOnMask & s < thresh & (1:numel(moveOnMask))' < peakIdx,1,'last');
            end
    end
    
    if(isempty(moveOnIdx))
        moveOnIdx = NaN;
    end
    
    
end

function [opts] = configureOpts(optsInput)

    opts = [];
    
    opts.PRE_TIME = -0.5;
    opts.POST_TIME = 1.2;
    
    opts.START_VAR = 'tgtOnTime';
    opts.END_VAR = 'endTime';
    
    opts.MOVE_START_OFFSET = 0;
    opts.MOVE_END_OFFSET = 100000;
    
    opts.MIN_VAR = 1.5;
    opts.THRESH_MULT = 0.5;
    opts.METHOD = 'peakAcceleration';
    
    opts.BUMP_MAGS = [];
    opts.STIM_CODES = [];
    
    opts.MAX_REACTION_TIME = 0.4; % in s
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
    
    opts.TARGET_DIRECTION = []; % initialize, but should not be set

end