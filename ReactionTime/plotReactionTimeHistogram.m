function [outputData,plots] = plotReactionTimeHistogram(reachData,cueInfo,opts)

    %% configure opts
    opts = configureOpts(opts);
    outputData = [];
    plots = [];
    
    %% for each cue idx, plot a histogram of the reaction times

    %% get bin counts
    cueIdx = 1;
    bE = opts.MIN_BIN:opts.BIN_SIZE:opts.MAX_BIN;
    %  each cue
    for cue = 1:numel(reachData.reactionTimes)
        bC(cueIdx,:) = histcounts(reachData.reactionTimes{cue},bE)/numel(reachData.reactionTimes{cue});
        cueIdx = cueIdx + 1;
    end
    
    
    %% make plot
    plots = figure();
    hold on
    for cueIdx = 1:size(bC,1)
        plot(bE(1:end-1)+mode(diff(bE))/2,bC(cueIdx,:),'-','color',opts.COLORS{cueIdx},'linewidth',opts.LINE_WIDTH)
    end
    
    %% set output data
    
    
end


function [opts] = configureOpts(optsInput)

    opts = [];
    
    opts.MIN_BIN = 0.17;
    opts.MAX_BIN = 0.27;
    opts.BIN_SIZE = 0.001;
    
    opts.LINE_WIDTH = 1.5;
    
    opts.COLORS = {'r',[0 0.5 0],'b','k','m',[0.5,0.5,0.2]};

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