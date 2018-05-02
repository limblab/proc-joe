function [f] = plotReactionTimeMeans(reachData,opts)
    % this function is meant to be run after plotReactionTimeHistogram.
    % inputData is expected to be the outputData from
    % plotReactionTimeHistogram
    
    %% configure opts
    opts = configureOpts(opts);
    
    %% plot bump data against bump magnitude
    
    bumpIndices = find(~isnan(reachData.bumpMags) & ~cellfun(@isempty,reachData.reactionTimes)')';
    
    % fit with a decaying exponential
    f = [];
    [f.fitObj,f.gof] = fit(reachData.bumpMags(bumpIndices),cellfun(@mean,reachData.reactionTimes(bumpIndices))','a*exp(b*x)+c','startPoint',[0,0,0.15]);
    
    % plot decaying exponential
    figure();
    hold on
    xData = linspace(min(reachData.bumpMags(bumpIndices)*0.9),max(reachData.bumpMags(bumpIndices)*1.1),100);
    yData = f.fitObj.a*exp(f.fitObj.b*xData)+f.fitObj.c;
    plot(xData,yData,'k--','linewidth',opts.LINE_WIDTH);

    % plot data points
    for bM = bumpIndices
        plot(reachData.bumpMags(bM),mean(reachData.reactionTimes{bM}),'.','color','k','markersize',opts.MARKER_SIZE)
    end
    % plot stim data against stim parameters
    
    % plot all conditions against probability of detection
end

function [opts] = configureOpts(optsInput)

    opts = [];
    
    opts.LINE_WIDTH = 1.5;
    
    opts.COLORS = {'r',[0 0.5 0],'b','k','m',[0.5,0.5,0.2]};
    opts.MARKER_SIZE = 20;
    
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