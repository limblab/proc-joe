function [reachPlot] = plotReachesTD(td,opts)

    %% configure opts
    opts = configureOpts(opts);
    
    %% extract td idxs for the plot
    if(~isempty(opts.BUMP_MAGS) && ~isempty(opts.STIM_CODES))
        td = td(isEqual([td.bumpMagnitude],opts.BUMP_MAGS) | isEqual([td.stimCode],opts.STIM_CODES));
    elseif(opts.KEEP_ONLY_VISUAL_TRIALS)
        td = td([td.isVisualTrial] == 1);
    elseif(~isempty(opts.BUMP_MAGS))
        td = td(isEqual([td.bumpMagnitude],opts.BUMP_MAGS));
    elseif(~isempty(opts.STIM_CODES))
        td = td(isEqual([td.stimCode],opts.STIM_CODES));
    end
    
    td = td(~isnan([td.idx_movement_on]));
    if(opts.RANDOM)
        td_reachPlot = datasample(td,min(numel(td),opts.MAX_PLOT),'replace',false);
    else
        td_reachPlot = td(3:3-1+min(numel(td),opts.MAX_PLOT));
    end
%     td_reachPlot = td(1:opts.MAX_PLOT);
    %% plot reaches with movement onset label
    figure();
    for i = 1:numel(opts.WHICH_IDX)
        subplot(numel(opts.WHICH_IDX),1,i);
        for r = 1:numel(td_reachPlot)
            t = ((1:size(td_reachPlot(r).pos,1)))*td_reachPlot(r).bin_size - td_reachPlot(r).goCueTime + 2*td_reachPlot(r).bin_size; % handle off by 1 error
            s = sum(td_reachPlot(r).(opts.WHICH_FIELD)(:,:).*[cos(opts.DIR*pi/180),sin(opts.DIR*pi/180)],2)';
            
            plot(t,s,'linewidth',opts.LINE_WIDTH,'color',opts.COLOR)

            % plot a dot for the reaction time
            hold on
            plot(t(td_reachPlot(r).idx_movement_on),s(td_reachPlot(r).idx_movement_on),'r.','markersize',opts.MARKER_SIZE);
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
    
    opts.DIR = 0;
    
    opts.XLIM = [-0.2,0.5];
    opts.YLIM = [];
    
    opts.MAX_PLOT = 10;
    
    opts.ZERO_MARKER = 'goCueTime'; % goCueTime is only supported option currently
    opts.WHICH_FIELD = 'vel';
    opts.WHICH_IDX = 1;
    opts.LINE_WIDTH = 1.5;
    opts.MARKER_SIZE = 20;
    
    opts.BUMP_MAGS = [];
    opts.STIM_CODES = [];
    opts.KEEP_ONLY_VISUAL_TRIALS = 0;
    
    opts.RANDOM = 1;
    opts.COLOR = 'k';
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