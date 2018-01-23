function [figureHandles] = plotHeatmaps(arrayData,opts)

    %% configure opts and set default values
    opts = configureOpts(opts);

    

end



function [opts] = configureOpts(optsInput)

    opts.BASELINE_PRE_TIME = -10/1000;
    opts.BASELINE_POST_TIME = -3/1000;
    opts.STIM_PRE_TIME = 1.2/1000;
    opts.STIM_POST_TIME = 5/1000;
    
    opts.STIM_ELECTRODE_PLOT = 1;
    opts.WAVEFORM_TYPES_PLOT = 1;
    
    opts.MAX_RATIO = 8;
    opts.MIN_RATIO = 1;
    opts.LOG_SCALE = 0;
    
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = '';
    opts.FIGURE_PREFIX = '';
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