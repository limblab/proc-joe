function [reachPlot] = plotReaches(reachData,cueInfo,opts)

    %% configure opts
    opts = configureOpts(opts);
    
    %% find a set of reaches to plot based on opts.NUM_PLOT
    reachIdx = datasample(find(~isnan(reachData.reactionTime)),min(opts.NUM_PLOT,sum(~isnan(reachData.reactionTime))),'replace',false);

    %% for each idx in reachIdx, plot the reach and a marker at the reaction time point
    figure();
    hold on
    for r = reachIdx'
        t = reachData.kin(r).t;
        if(strcmpi(opts.ZERO_MARKER,'goCueTime')==1) % shift t over so that goCueTime is 0
            t = t - reachData.goCueTime(r);
        end
        plot(t,reachData.kin(r).(opts.PLOT_VAR),'k','linewidth',opts.LINE_WIDTH)
        
        % plot a dot for the reaction time
        rt = reachData.reactionTime(r);
        if(~strcmpi(opts.ZERO_MARKER,'goCueTime')==1)
            rt = rt + reachData.goCueTime(r);
        end
        dataIdx = find(t >= rt,1,'first');
        
        plot(rt,reachData.kin(r).(opts.PLOT_VAR)(dataIdx),'r.','markersize',opts.MARKER_SIZE);
    end
    
    xlim([opts.PRE_TIME,opts.POST_TIME])
end



function [opts] = configureOpts(optsInput)

    opts = [];
    
    opts.PRE_TIME = -0.2;
    opts.POST_TIME = 0.5;
    
    opts.NUM_PLOT = 10;
    
    opts.ZERO_MARKER = 'goCueTime'; % goCueTime and an else are the only supported options currently
    opts.PLOT_VAR = 'vx';
    opts.LINE_WIDTH = 1.5;
    opts.MARKER_SIZE = 20;
    
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