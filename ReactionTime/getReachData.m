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
            end
            % find reaction time
            reachData.reactionTime(tr) = reachData.kin(tr).moveOnTime - reachData.goCueTime(tr);
            
        end
    end
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
    moveOnMask(1:find(kin.t > cueTime,1,'first')+opts.MOVE_START_OFFSET) = 0;
    
    % with remaining data, find peak vel baesd on accel crossing zero,
    % finds first peak
    s = kin.vx;
    ds = [0;diff(s)];
    dds = [0;diff(ds)];
    peaks = [dds(1:end-1) > 0 & dds(2:end) < 0; 0];
    peakIdx = find(peaks & ds > opts.MIN_DS & moveOnMask,1,'first');
    
    if(~isempty(peakIdx))
        % find latest time at which the acceleration is half maximum
        thresh = kin.ax(peakIdx)/2;
        moveOnMask = moveOnMask & kin.ax > thresh;
        
        moveOnIdx = find(moveOnMask,1,'first');
    end
    
    
end

function [opts] = configureOpts(optsInput)

    opts = [];
    
    opts.PRE_TIME = -1.2;
    opts.POST_TIME = 1.2;
    
    opts.START_VAR = 'tgtOnTime';
    opts.END_VAR = 'endTime';
    
    opts.MOVE_START_OFFSET = 0;
    opts.MIN_DS = 0.2;
    
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