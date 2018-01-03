function [ figureHandles ] = plotWaveformsStim( cds, NEURON_NUMBER, opts)
%%% need to do save figures, titles, and other plot format stuff
%%% also align waves?


    %% configure opts and set default values
    opts = configureOpts(opts);

    %% get number of chans, chan list, and waveform list
    if(any(isfield(cds,'waveforms')))
        NUM_WAVEFORM_TYPES = numel(unique(cds.waveforms.waveSent));
    else
        NUM_WAVEFORM_TYPES = 1;
    end

    if(any(isfield(cds.waveforms,'chanSent')))
        NUM_CHANS = numel(unique(cds.waveforms.chanSent));
        CHAN_LIST = unique(cds.waveforms.chanSent);
    else
        CHAN_LIST = opts.STIM_ELECTRODE;
        NUM_CHANS = 1;
    end
    figureHandles = {};
    
    %% loop over conditions and plot waveforms
    for chan = 1:NUM_CHANS
        for wave = 1:NUM_WAVEFORM_TYPES
            % get waveform indexes for plotting
            spikeMask.nearArtifact = zeros(numel(cds.units(NEURON_NUMBER).spikes.ts,1));
            spikeMask.awayArtifact = zeros(numel(cds.units(NEURON_NUMBER).spikes.ts,1));
            
            for st = 1:numel(cds.stimOn)-1 % not important to do all of them
                if(cds.waveforms.chanSent(st) == CHAN_LIST(chan) && cds.waveforms.waveSent(st) == wave)
                    spikeMask.nearArtifact = spikeMask.nearArtifact | (cds.units(NEURON_NUMBER).spikes.ts > cds.stimOn(st) & ...
                        cds.units(NEURON_NUMBER).spikes.ts < cds.stimOn(st) + opts.TIME_AFTER_STIMULATION_ARTIFACT);

                    spikeMask.awayArtifact = spikeMask.awayArtifact | (cds.units(NEURON_NUMBER).spikes.ts > cds.stimOn(st) + opts.TIME_AFTER_STIMULATION_NO_ARTIFACT & ...
                        cds.units(NEURON_NUMBER).spikes.ts < cds.stimOn(st+1));
                end
            end
            
            waveIdx.nearArtifact = find(spikeMask.nearArtifact);
            waveIdx.awayArtifact = find(spikeMask.awayArtifact);
            
            if(opts.RANDOM_SAMPLE) % rabdomly sample
                waveIdx.nearArtifact = datasample(waveIdx.nearArtifact,min(opts.MAX_WAVES_PLOT,numel(waveIdx.nearArtifact)),'replace',false);
                waveIdx.awayArtifact = datasample(waveIdx.awayArtifact,min(opts.MAX_WAVES_PLOT,numel(waveIdx.awayArtifact)),'replace',false);
            else % grab first set
                waveIdx.nearArtifact = waveIdx.nearArtifact(1:min(opts.MAX_WAVES_PLOT,numel(waveIdx.nearArtifact)));
                waveIdx.awayArtifact = waveIdx.awayArtifact(1:min(opts.MAX_WAVES_PLOT,numel(waveIdx.awayArtifact)));
            end
            
            % get rawIdxs
            rawIdx.nearArtifact = getRawDataIdx(cds.units(NEURON_NUMBER).spikes.ts(waveIdx.nearArtifact),zeros(size(waveIdx.nearArtifact)) + double(cds.units(NEURON_NUMBER).chan),cds.rawData.ts,cds.rawData.elec);
            rawIdx.awayArtifact = getRawDataIdx(cds.units(NEURON_NUMBER).spikes.ts(waveIdx.awayArtifact),zeros(size(waveIdx.awayArtifact)) + double(cds.units(NEURON_NUMBER).chan),cds.rawData.ts,cds.rawData.elec);
            
            rawIdx.nearArtifact(rawIdx.nearArtifact == -1) = [];
            rawIdx.awayArtifact(rawIdx.awayArtifact == -1) = [];
            
            % we have all the indexes we need, now plot the waveforms
            if(sum(opts.PLOT_FILTERED == 1) > 0) % plot filtered waveforms (cds.units...spikes)
                figureHandles{end+1} = figure(); % filtered near artifact
                xDataFiltered = ((1:size(cds.units(NEURON_NUMBER).spikes{:,2:end},2))-1)/30;
                
                plot(xDataFiltered,cds.units(NEURON_NUMBER).spikes.wave(waveIdx.nearArtifact,:)-mean(cds.units(NEURON_NUMBER).spikes.wave(waveIdx.nearArtifact,end-10:end),2))
                ylim([opts.YLIM_FILTERED])
                
                figureHandles{end+1} = figure(); % filtered away from artifact
                plot(xDataFiltered,cds.units(NEURON_NUMBER).spikes.wave(waveIdx.awayArtifact,:)-mean(cds.units(NEURON_NUMBER).spikes.wave(waveIdx.awayArtifact,end-10:end),2))
                ylim([opts.YLIM_FILTERED])
            end
            if(sum(opts.PLOT_FILTERED == 0) > 0) % plot raw waveforms (cds.rawData)
                figureHandles{end+1} = figure(); % raw near artifact
                xDataRaw = ((1:size(cds.rawData.waveforms,2))-1)/30;
                plot(xDataRaw,cds.rawData.waveforms(rawIdx.nearArtifact,:) - mean(cds.rawData.waveforms(rawIdx.nearArtifact,end-10:end),2))
                
                figureHandles{end+1} = figure(); % filtered away from artifact
                plot(xDataRaw,cds.rawData.waveforms(rawIdx.awayArtifact,:) - mean(cds.rawData.waveforms(rawIdx.awayArtifact,end-10:end),2))
            end
            
            
            
        end
    end
    
    
end % end function


function [opts] = configureOpts(optsInput)

    opts.MAX_WAVES_PLOT = 100;
    opts.RANDOM_SAMPLE = 1;
    
    opts.TIME_AFTER_STIMULATION_ARTIFACT = 5/1000;
    opts.TIME_AFTER_STIMULATION_NO_ARTIFACT = 20/1000;
    
    opts.PLOT_FILTERED = [0,1];
    
    opts.PLOT_TITLE = 1;
    opts.TITLE_TO_PLOT = '';
    
    opts.ALIGN_WAVES = 1;
    opts.ADJUST_AMPLITUDE = 1;
    
    opts.NEURON_NUMBER = 1;
    
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
end % end function