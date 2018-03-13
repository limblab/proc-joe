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
            
            if(opts.RANDOM_SAMPLE) % randomly sample
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
                
                waveforms.nearArtifact = cds.units(NEURON_NUMBER).spikes.wave(waveIdx.nearArtifact,:);
                waveforms.awayArtifact = cds.units(NEURON_NUMBER).spikes.wave(waveIdx.awayArtifact,:);
                
                if(opts.ALIGN_WAVES)
                    waveforms.nearArtifact = alignWaves(waveforms.nearArtifact,'max','filtered',opts);
                    waveforms.awayArtifact = alignWaves(waveforms.awayArtifact,'max','filtered',opts);
                end
                                
                opts.FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_nn',num2str(NEURON_NUMBER),'_chan',num2str(cds.units(NEURON_NUMBER).chan),'stimChan',num2str(CHAN_LIST(chan)),'_waveNum',num2str(wave),'wavesFilteredNearArtifact');
                figureHandles{end+1} = plotWaveforms(waveforms.nearArtifact,opts);
                
                opts.FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_nn',num2str(NEURON_NUMBER),'_chan',num2str(cds.units(NEURON_NUMBER).chan),'stimChan',num2str(CHAN_LIST(chan)),'_waveNum',num2str(wave),'wavesFilteredAwayArtifact');
                figureHandles{end+1} = plotWaveforms(waveforms.awayArtifact,opts);
                
            end
            if(sum(opts.PLOT_FILTERED == 0) > 0) % plot raw waveforms (cds.rawData)
                
                
                waveforms.nearArtifact = cds.rawData.waveforms(rawIdx.nearArtifact,opts.PRE_WAVE_RAW:opts.PRE_WAVE_RAW + opts.WAVE_LENGTH_RAW);
                waveforms.awayArtifact = cds.rawData.waveforms(rawIdx.awayArtifact,opts.PRE_WAVE_RAW:opts.PRE_WAVE_RAW + opts.WAVE_LENGTH_RAW);
                
                if(opts.ALIGN_WAVES)
                    waveforms.nearArtifact = alignWaves(waveforms.nearArtifact,'min','raw',opts);
                    waveforms.awayArtifact = alignWaves(waveforms.awayArtifact,'min','raw',opts);
                end
                                
                opts.FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_nn',num2str(NEURON_NUMBER),'_chan',num2str(cds.units(NEURON_NUMBER).chan),'stimChan',num2str(CHAN_LIST(chan)),'_waveNum',num2str(wave),'wavesRawNearArtifact');
                figureHandles{end+1} = plotWaveforms(waveforms.nearArtifact,opts);
                
                opts.FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_nn',num2str(NEURON_NUMBER),'_chan',num2str(cds.units(NEURON_NUMBER).chan),'stimChan',num2str(CHAN_LIST(chan)),'_waveNum',num2str(wave),'wavesRawAwayArtifact');
                figureHandles{end+1} = plotWaveforms(waveforms.awayArtifact,opts);
                
            end
            
            
            
        end
    end
    
    
end % end function






function [waveformsAligned] = alignWaves(waveforms,maxOrMin,rawOrFiltered,opts)

    waveformsAligned = waveforms;
    
    alignOffset = opts.ALIGN_OFFSET;
    middleIdx = opts.MIDDLE_IDX;
    
    if(strcmpi(rawOrFiltered,'raw') == 1)
        middleIdx = opts.MIDDLE_IDX_RAW;
        if(opts.ADJUST_FOR_BIT_ERROR)
            waveforms = waveforms/0.254;
        end
    end
    
    % find max or min idx
    if(strcmpi(maxOrMin,'min') == 1)
        [~,peakIdx] = min(waveforms(:,alignOffset:end),[],2);
    elseif(strcmpi(maxOrMin,'max') == 1)
        [~,peakIdx] = max(waveforms(:,alignOffset:end),[],2);
    else
        error('max or min not specified, nothing done')
    end
    
    peakIdx = peakIdx + alignOffset;

    % shift waves accordingly
    shiftIdx = peakIdx - middleIdx;

    for w = 1:size(waveforms,1)
        waveformsAligned(w,:) = circshift(waveforms(w,:),-1*shiftIdx(w));
        waveformsAligned(w,end-shiftIdx(w):end) = 0;
    end
                    
end

function [figureHandle] = plotWaveforms(waveforms,opts)
    figureHandle = figure();  
    
    % deal with colors
    if(~isempty(opts.COLOR_RANGE))
        colorOrder = [linspace(opts.COLOR_RANGE(1,1),opts.COLOR_RANGE(2,1),size(waveforms,1))',...
            linspace(opts.COLOR_RANGE(1,2),opts.COLOR_RANGE(2,2),size(waveforms,1))',...
            linspace(opts.COLOR_RANGE(1,3),opts.COLOR_RANGE(2,3),size(waveforms,1))'];
              
        set(gca,'colorOrder',colorOrder);
        hold all
    end
    
    % plot
    if(~isempty(waveforms))
        xData = ((1:size(waveforms,2))-1)/30;
        plot(xData,waveforms-mean(waveforms(:,end-10:end),2))
        ylim(opts.YLIM)

        formatForLee(gcf)
    end
    
%     if(opts.PLOT_TITLE)
%         title(opts.TITLE_TO_PLOT)
%     end
    
    % deal with saving plot
    if(opts.FIGURE_SAVE && strcmpi(opts.FIGURE_NAME,'')~=1 && strcmpi(opts.FIGURE_DIR,'')~=1)
        saveFiguresLIB(gcf,opts.FIGURE_DIR,strcat(opts.FIGURE_NAME));
    end
    
    
end



function [opts] = configureOpts(optsInput)

    opts.ADJUST_FOR_BIT_ERROR = 1;

    opts.MAX_WAVES_PLOT = 100;
    opts.RANDOM_SAMPLE = 1;
    
    opts.TIME_AFTER_STIMULATION_ARTIFACT = 5/1000;
    opts.TIME_AFTER_STIMULATION_NO_ARTIFACT = 20/1000;
    
    opts.PLOT_FILTERED = [0,1];
    
    opts.PLOT_TITLE = 0; % not implemented
    opts.TITLE_TO_PLOT = '';
    
    opts.ALIGN_WAVES = 1;
    opts.PRE_WAVE = 22;
    opts.POST_WAVE = 25;
    opts.PRE_WAVE_RAW = 28;
    opts.WAVE_LENGTH_RAW = 48;
    
    opts.ALIGN_OFFSET = 18;
        
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = '';
    opts.FIGURE_PREFIX = '';
    
    opts.YLIM = [-500,500];
    opts.COLOR_RANGE = [0.0,0.0,0.0;0.5,0.5,0.5];
    
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
    opts.MIDDLE_IDX = opts.PRE_WAVE + 1;
    opts.MIDDLE_IDX_RAW = opts.MIDDLE_IDX + 4;
    
end % end function