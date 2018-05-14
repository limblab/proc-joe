function [reachPlot] = plotReachesTD(td,opts)

    %% configure opts
    opts = configureOpts(opts);
    
    %% extract td idxs for the plot
    if(~isempty(opts.BUMP_MAGS) && ~isempty(opts.STIM_CODES))
        td = td(isEqual([td.bumpMagnitude],opts.BUMP_MAGS) | isEqual([td.stimCode],opts.STIM_CODES));
    elseif(~isempty(opts.BUMP_MAGS))
        td = td(isEqual([td.bumpMagnitude],opts.BUMP_MAGS));
    elseif(~isempty(opts.STIM_CODES))
        td = td(isEqual([td.stimCode],opts.STIM_CODES));
    end
    td_reachPlot = datasample(td,min(numel(td),opts.MAX_PLOT),'replace',false);
    
    %% plot reaches with movement onset label
    figure();
    for i = 1:numel(opts.WHICH_IDX)
        subplot(numel(opts.WHICH_IDX),1,i);
        for r = 1:numel(td_reachPlot)
            t = ((1:size(td_reachPlot(r).pos,1))-td_reachPlot(r).idx_goCueTime)*td_reachPlot(r).bin_size;

            plot(t,td_reachPlot(r).(opts.WHICH_FIELD)(:,opts.WHICH_IDX(i)),'k','linewidth',opts.LINE_WIDTH)

            % plot a dot for the reaction time
            hold on
            plot(t(td_reachPlot(r).idx_movement_on),td_reachPlot(r).(opts.WHICH_FIELD)(td_reachPlot(r).idx_movement_on,opts.WHICH_IDX(i)),'r.','markersize',opts.MARKER_SIZE);
        end
        if(~isempty(opts.XLIM))
            xlim(opts.XLIM);
        end
        if(~isempty(opts.YLIM))
            ylim(opts.YLIM);
        end
        formatForLee(gcf);
    end
    

end


function [opts] = configureOpts(optsInput)

    opts = [];
    
    opts.XLIM = [-0.2,0.5];
    opts.YLIM = [];
    
    opts.MAX_PLOT = 10;
    
    opts.ZERO_MARKER = 'goCueTime'; % goCueTime and an else are the only supported options currently
    opts.WHICH_FIELD = 'vel';
    opts.WHICH_IDX = 1;
    opts.LINE_WIDTH = 1.5;
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
    

end