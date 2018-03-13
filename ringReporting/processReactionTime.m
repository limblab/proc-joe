function [output_data] = processReactionTime(cds,opts)

    %% configure opts
    opts = configureOpts(opts);
    figureHandles = {};

    %% get reaction time
    
    reactionTimes = getReactionTime(cds,opts);

    output_data = 0;
end


function [rt] = getReactionTime(cds,opts)
    rt = 0;
    figure
    % loop through each trial
    for tr = 1:30%size(cds.trials,1)
        if(cds.trials.result(tr) == 'R' || cds.trials.result(tr) == 'F') % if trial was attempted
            kinIdx = [find(cds.trials.goCueTime(tr) < cds.kin.t,1,'first'), find(cds.trials.endTime(tr) < cds.kin.t,1,'first')];
            
            plot(cds.kin.t(kinIdx(1):kinIdx(2)) - cds.kin.t(kinIdx(1)),sqrt(cds.kin.vx(kinIdx(1):kinIdx(2)).^2 + cds.kin.vy(kinIdx(1):kinIdx(2)).^2));
            hold on
            
        end
    end

end


function [opts] = configureOpts(optsInput)

    opts = [];
    opts.MAKE_FIGURES = 1;
    opts.FONT_SIZE = 14;
    opts.LINE_WIDTH = 2;
    opts.CENTER_X = 3;
    opts.CENTER_Y = -33;
    
    
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_PREFIX = '';
    opts.FIGURE_DIR = '';
    
    %% check if in opts and optsInput, overwrite if so
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