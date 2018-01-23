function [ figureHandle ] = plotInterspikeIntervalHistogram( cds, NEURON_NUMBER, opts )

    %% configure opts and set default values
    opts = configureOpts(opts);
    
    figureHandle = figure();
    %% get data
    interSpikeInterval = diff(cds.units(NEURON_NUMBER).spikes.ts);
    interSpikeInterval = interSpikeInterval*1000; % seconds to ms

    %% bin and plot data
    binEdges = 0:opts.BIN_SIZE*1000:max(opts.XLIM);
    [binCounts,~] = histcounts(interSpikeInterval,binEdges);
    bar(binEdges(1:end-1)+(binEdges(2)-binEdges(1))/2,binCounts)
    ylabel('Spike count')
    xlabel('Time (ms)')
    formatForLee(gcf)
    %% display things (namely percentage below 1.7ms)
    if(opts.DISPLAY_TEXT)
        percentageBelow = sum(interSpikeInterval < 1.7)/numel(interSpikeInterval)*100;
        disp(percentageBelow)
    end

    if(opts.FIGURE_SAVE && strcmpi(opts.FIGURE_NAME,'')~=1 && strcmpi(opts.FIGURE_DIR,'')~=1)
        saveFiguresLIB(figHandle,opts.FIGURE_DIR,strcat(opts.FIGURE_PREFIX,'_nn',num2str(NEURON_NUMBER),'_chan',num2str(cds.units(NEURON_NUMBER).chan),'_ISI'));
    end

end


function [opts] = configureOpts(optsInput)

    opts.XLIM = [0,20];
    opts.BIN_SIZE = 1/1000;
    opts.DISPLAY_TEXT = 1;
    
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
    
    %% constants that depend on variables above
    
end % end function