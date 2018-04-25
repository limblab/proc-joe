function [behaviorData] = getReactionTimeSuccessRates(cds,cueInfo,opts)

    % this function will get the 4 things:
    % no cue trials : % stay still and % false start
    % cue trials: % move, % no move based split up by cue
    
    
    %% configure opts
    opts = configureOpts(opts);
    
    
    %% find no cue trials
    noCueIdx = find(~isnan([cds.trials.tgtOnTime]) & isnan([cds.trials.bumpTime]));
    
    behaviorData.falseStartRate = sum([cds.trials.result(noCueIdx)] == 'A')/numel(noCueIdx);
    behaviorData.stayStillRate = sum([cds.trials.result(noCueIdx)] == 'R')/numel(noCueIdx);
    
    %% find success rate for trials with a cue based on cue type
    c = 1;
    bumpMags = unique(cueInfo.bumpMag(~isnan(cueInfo.bumpMag)));
    stimCodes = unique(cueInfo.stimCode(~isnan(cueInfo.stimCode)));
    %  each bump magnitude
    for bM = 1:numel(bumpMags)
        cueIdx = find(cueInfo.bumpMag == bumpMags(bM));
        behaviorData.reachSuccess(c) = sum(cds.trials.result(cueIdx) == 'R')/numel(cueIdx);
        behaviorData.reachNums(c) = numel(cueIdx);
        behaviorData.bumpMag(c) = bumpMags(bM);
        behaviorData.stimCode(c) = NaN;
        
        c = c + 1;
    end
    
    %  each stim code
    for sC = 1:numel(stimCodes)
        behaviorData.reachSuccess(c) = sum(cds.trials.result(cueIdx) == 'R')/numel(cueIdx);
        behaviorData.reachNums(c) = numel(cueIdx);
        behaviorData.bumpMag(c) = NaN;
        behaviorData.stimCode(c) = stimCodes(sc);
        c = c + 1;
    end


end


function [opts] = configureOpts(optsInput)

    opts = [];
    
    
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