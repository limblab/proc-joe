function [ figureHandle ] = plotInterspikeIntervalHistogram( cds, NEURON_NUMBER, opts )

    %% configure opts and set default values
    opts = configureOpts(opts);
    
    figureHandle = figure();
    %% get data
    spikeTimes = cds.units(NEURON_NUMBER).spikes.ts;
    
    % clean out spike times
    spikeTimesKeep = spikeTimes(spikeTimes<10);
    
    for st = 1:numel(cds.stimOn)-1
        spikeTimesKeep = [spikeTimesKeep; spikeTimes(spikeTimes > cds.stimOn(st) + 0.08 & spikeTimes < cds.stimOn(st+1))];
%         spikeTimesKeep = [spikeTimesKeep; spikeTimes(cds.waveforms.waveSent(st) == 2 & spikeTimes > cds.stimOn(st) & spikeTimes < cds.stimOn(st)+0.15)];
    end   
    spikeTimesKeep = [spikeTimesKeep; spikeTimes(spikeTimes > cds.stimOn(end) + 0.1)];
    spikeTimes = spikeTimesKeep;

    interSpikeInterval = diff(spikeTimes);
    interSpikeInterval = interSpikeInterval*1000; % seconds to ms
%     interSpikeInterval = interSpikeInterval(interSpikeInterval > 1.5);
    %% bin and plot data
    binEdges = 0:opts.BIN_SIZE*1000:max(opts.XLIM);
    [binCounts,~] = histcounts(interSpikeInterval,binEdges);
    bar(binEdges(1:end-1)+(binEdges(2)-binEdges(1))/2,binCounts)
    ylabel('Spike count')
    xlabel('Time (ms)')
    formatForLee(gcf)
    %% display things (namely percentage below 1.7ms)
    if(opts.DISPLAY_TEXT)
        percentageBelow = sum(interSpikeInterval < 1.2)/numel(interSpikeInterval)*100;
        disp(percentageBelow)
    end

    if(opts.FIGURE_SAVE && strcmpi(opts.FIGURE_PREFIX,'')~=1 && strcmpi(opts.FIGURE_DIR,'')~=1)
        saveFiguresLIB(gcf,opts.FIGURE_DIR,strcat(opts.FIGURE_PREFIX,'_nn',num2str(NEURON_NUMBER),'_chan',num2str(cds.units(NEURON_NUMBER).chan),'_ISI'));
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