function [outputData,plots] = plotReactionTimeHistogram(reachData,cueInfo,opts)

    %% configure opts
    opts = configureOpts(opts);
    outputData = [];
    plots = [];
    
    %% for each cue idx, plot a histogram of the reaction times
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

    %% get bin counts
    cueIdx = 1;
    bE = opts.MIN_BIN:opts.BIN_SIZE:opts.MAX_BIN;
    totalCount = 0;
    %  each bump magnitude
    for bM = 1:numel(bumpMags)
        reactionTimes = reachData.reactionTime(find(isEqual(cueInfo.bumpMag,bumpMags(bM)) & ~isnan(reachData.reactionTime)));
        bC(cueIdx,:) = histcounts(reactionTimes,bE)/sum(isEqual(cueInfo.bumpMag,bumpMags(bM))&~isnan(reachData.reactionTime));
        cueIdx = cueIdx + 1;
        totalCount = totalCount + sum(isEqual(cueInfo.bumpMag,bumpMags(bM))&~isnan(reachData.reactionTime));
    end
    
    %  each stim code
    for sC = 1:numel(stimCodes)
        reactionTimes = reachData.reactionTime(isEqual(cueInfo.stimCode,stimCodes(sC)&~isnan(reachData.reactionTime)));
        bC(cueIdx,:) = histcounts(reactionTimes,bE)/sum(isEqual(cueInfo.stimCode,stimCodes(sC))&~isnan(reachData.reactionTime));
        cueIdx = cueIdx + 1;
    end
    
    %% make plot
    plots = figure();
    hold on
    for cueIdx = 1:size(bC,1)
        plot(bE(1:end-1)+mode(diff(bE))/2,bC(cueIdx,:),'-','color',opts.COLORS{cueIdx},'linewidth',opts.LINE_WIDTH)
    end
    
    %% set output data
    outputData.bC = bC;
    outputData.bE = bE;
    outputData.bumpMags = [bumpMags;NaN(numel(stimCodes),1)];
    outputData.stimCodes = [NaN(numel(bumpMags),1);stimCodes];
    
end

function [out] = isEqual(data1,data2,thresh)

    out = data1 < data2 + eps & data1 > data2 - eps;

end

function [opts] = configureOpts(optsInput)

    opts = [];
    
    opts.MIN_BIN = 0.17;
    opts.MAX_BIN = 0.27;
    opts.BIN_SIZE = 0.001;
    
    opts.LINE_WIDTH = 1.5;
    
    opts.COLORS = {'r',[0 0.5 0],'b','k','m',[0.5,0.5,0.2]};
    
    opts.BUMP_MAGS = [];
    opts.STIM_CODES = [];
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
    
    % check for row vectors instead of column vectors
    if(size(opts.BUMP_MAGS,2) > 1)
        opts.BUMP_MAGS = opts.BUMP_MAGS';
    end
    if(size(opts.STIM_CODES,2) > 1)
        opts.STIM_CODES = opts.STIM_CODES';
    end

end