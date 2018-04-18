function [ figureHandle, outputData ] = plotInterspikeIntervalHistogram( cds,stimInfo NEURON_NUMBER, opts )

    %% configure opts and set default values
    opts = configureOpts(opts);
    
    %% get data
    spikeTimes = cds.units(NEURON_NUMBER).spikes.ts;
    
    % clean out spike times
    ts_postStim = [];
    ts_noStim = [];
    
    for st = 1:numel(stimInfo.stimOn)-1
        ts_noStim = [ts_noStim; spikeTimes(spikeTimes > stimInfo.stimOn(st) + 0.05 & spikeTimes < stimInfo.stimOn(st+1))];
        ts_postStim = [ts_postStim; spikeTimes(spikeTimes > stimInfo.stimOn(st) & spikeTimes < stimInfo.stimOn(st)+0.02)];
    end   
    ts_noStim = [spikeTimes(spikeTimes < stimInfo.stimOn(1)); ts_noStim; spikeTimes(spikeTimes > stimInfo.stimOn(end) + 0.1)];

    interSpikeInterval{1} = diff(ts_noStim)*1000;
    interSpikeInterval{2} = diff(ts_postStim)*1000; % seconds to ms

    %% bin and plot data
    binEdges = 0:opts.BIN_SIZE*1000:max(max(interSpikeInterval{1}),max(interSpikeInterval{2}));
    for i = 1:numel(interSpikeInterval)
        figureHandle{i} = figure();

        [binCounts,~] = histcounts(interSpikeInterval{i},binEdges);
        bar(binEdges(1:end-1)+(binEdges(2)-binEdges(1))/2,binCounts)
        ylabel('Spike count')
        xlabel('Time (ms)')
        xlim([opts.XLIM])
        formatForLee(gcf)
        %% display things (namely percentage below 1.7ms)
        if(opts.DISPLAY_TEXT)
            percentageBelow = sum(interSpikeInterval{i} < 1.2)/numel(interSpikeInterval{i})*100;
            disp(percentageBelow)
        end

        if(opts.FIGURE_SAVE && strcmpi(opts.FIGURE_PREFIX,'')~=1 && strcmpi(opts.FIGURE_DIR,'')~=1)
            saveFiguresLIB(gcf,opts.FIGURE_DIR,strcat(opts.FIGURE_PREFIX,'_nn',num2str(NEURON_NUMBER),'_chan',num2str(cds.units(NEURON_NUMBER).chan),'_ISI'));
        end
    end

    
    outputData.interSpikeInterval = interSpikeInterval;
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