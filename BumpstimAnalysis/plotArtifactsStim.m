function [figureHandles] = plotArtifactsStim(cds,stimInfo,artifactData,NEURON_NUMBER,opts)

    %% configure opts and set default values
    opts = configureOpts(opts);

    %% get number of chans, chan list, and waveform list
    if(any(isfield(stimInfo,'waveforms')))
        NUM_WAVEFORM_TYPES = numel(unique(stimInfo.waveforms.waveSent));
    else
        NUM_WAVEFORM_TYPES = 1;
    end

    if(any(isfield(stimInfo.waveforms,'chanSent')))
        NUM_CHANS = numel(unique(stimInfo.waveforms.chanSent));
        CHAN_LIST = unique(stimInfo.waveforms.chanSent);
    else
        CHAN_LIST = opts.STIM_ELECTRODE;
        NUM_CHANS = 1;
    end
    figureHandles = {};
    if(opts.ARTIFACT_MULTIPLIER == -1)
        opts.ARTIFACT_MULTIPLIER = numel(stimInfo.stimOn)/size(artifactData.artifact,1); % no longer keep all artifacts, use this to shift indexes
    end
    %% loop over conditions
    for chan = 1:NUM_CHANS
        for wave = 1:NUM_WAVEFORM_TYPES
            if((isempty(opts.WAVEFORMS_TO_PLOT) || sum(wave==opts.WAVEFORMS_TO_PLOT) > 0) && (isempty(opts.CHANS_TO_PLOT) || sum(chan==opts.CHANS_TO_PLOT) > 0))
                % get artifact indices
                artifactIdx.nearArtifact = [];
                artifactIdx.awayArtifact = [];

                for st = 1:opts.ARTIFACT_MULTIPLIER:numel(stimInfo.stimOn)
                    if(stimInfo.waveforms.chanSent(st) == CHAN_LIST(chan) && stimInfo.waveforms.waveSent(st) == wave)
                        artifactMask = cds.units(NEURON_NUMBER).spikes.ts > artifactData.t((st-1)/opts.ARTIFACT_MULTIPLIER + 1) + 4/1000 & ...
                            cds.units(NEURON_NUMBER).spikes.ts < artifactData.t((st-1)/opts.ARTIFACT_MULTIPLIER + 1) + opts.TIME_AFTER_STIMULATION_ARTIFACT & ...
                            ~(cds.units(NEURON_NUMBER).spikes.ts > artifactData.t((st-1)/opts.ARTIFACT_MULTIPLIER + 1) & ...
                            cds.units(NEURON_NUMBER).spikes.ts < artifactData.t((st-1)/opts.ARTIFACT_MULTIPLIER + 1) +4/1000);
                        if(sum(artifactMask)==0) % no spikes present
                            artifactIdx.awayArtifact(end+1,1) = (st-1)/opts.ARTIFACT_MULTIPLIER + 1;
                        elseif(sum(artifactMask)>0) % spike present
                            artifactIdx.nearArtifact(end+1,1) = (st-1)/opts.ARTIFACT_MULTIPLIER + 1;
                        end
                    end
                end

                if(opts.RANDOM_SAMPLE) % randomly sample
                    artifactIdx.nearArtifact = datasample(artifactIdx.nearArtifact,min(opts.MAX_WAVES_PLOT*opts.ROW_SUBPLOT*opts.COL_SUBPLOT,numel(artifactIdx.nearArtifact)),'replace',false);
                    artifactIdx.awayArtifact = datasample(artifactIdx.awayArtifact,min(opts.MAX_WAVES_PLOT*opts.ROW_SUBPLOT*opts.COL_SUBPLOT,numel(artifactIdx.awayArtifact)),'replace',false);
                else % grab first set
                    artifactIdx.nearArtifact = artifactIdx.nearArtifact(1:min(opts.MAX_WAVES_PLOT*opts.ROW_SUBPLOT*opts.COL_SUBPLOT,numel(artifactIdx.nearArtifact)));
                    artifactIdx.awayArtifact = artifactIdx.awayArtifact(1:min(opts.MAX_WAVES_PLOT*opts.ROW_SUBPLOT*opts.COL_SUBPLOT,numel(artifactIdx.awayArtifact)));
                end

                % we have all the indexes we need, now plot the waveforms
                if(sum(opts.PLOT_FILTERED == 1) > 0) % plot filtered waveforms (cds.units...spikes)

                    waveforms.nearArtifact = squeeze(artifactData.artifact(artifactIdx.nearArtifact,cds.units(NEURON_NUMBER).chan,:));
                    waveforms.awayArtifact = squeeze(artifactData.artifact(artifactIdx.awayArtifact,cds.units(NEURON_NUMBER).chan,:));

                    waveforms.nearArtifact = filterArtifactData(waveforms.nearArtifact,opts);
                    waveforms.awayArtifact = filterArtifactData(waveforms.awayArtifact,opts);

                    opts.FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_nn',num2str(NEURON_NUMBER),'_chan',num2str(cds.units(NEURON_NUMBER).chan),'stimChan',num2str(CHAN_LIST(chan)),'_waveNum',num2str(wave),'artifactsFilteredNear');
                    opts.TITLE_TO_PLOT = strcat('Filtered-Near, stimChan ',num2str(CHAN_LIST(chan)),', wave ', num2str(wave));
                    figureHandles{end+1} = plotArtifacts(waveforms.nearArtifact,opts);

                    opts.FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_nn',num2str(NEURON_NUMBER),'_chan',num2str(cds.units(NEURON_NUMBER).chan),'stimChan',num2str(CHAN_LIST(chan)),'_waveNum',num2str(wave),'artifactsFilteredAway');
                    opts.TITLE_TO_PLOT = strcat('Filtered-Away, stimChan ',num2str(CHAN_LIST(chan)),', wave ', num2str(wave));
                    figureHandles{end+1} = plotArtifacts(waveforms.awayArtifact,opts);

                end
                if(sum(opts.PLOT_FILTERED == 0) > 0) % plot raw waveforms (cds.rawData)

                    waveforms.nearArtifact = squeeze(artifactData.artifact(artifactIdx.nearArtifact,cds.units(NEURON_NUMBER).chan,:));
                    waveforms.awayArtifact = squeeze(artifactData.artifact(artifactIdx.awayArtifact,cds.units(NEURON_NUMBER).chan,:));

                    opts.FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_nn',num2str(NEURON_NUMBER),'_chan',num2str(cds.units(NEURON_NUMBER).chan),'stimChan',num2str(CHAN_LIST(chan)),'_waveNum',num2str(wave),'artifactsRawNear');
                    opts.TITLE_TO_PLOT = strcat('Raw-Near, stimChan ',num2str(CHAN_LIST(chan)),', wave ', num2str(wave));
                    figureHandles{end+1} = plotArtifacts(waveforms.nearArtifact,opts);

                    opts.FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_nn',num2str(NEURON_NUMBER),'_chan',num2str(cds.units(NEURON_NUMBER).chan),'stimChan',num2str(CHAN_LIST(chan)),'_waveNum',num2str(wave),'artifactsRawAway');
                    opts.TITLE_TO_PLOT = strcat('Raw-Away, stimChan ',num2str(CHAN_LIST(chan)),', wave ', num2str(wave));
                    figureHandles{end+1} = plotArtifacts(waveforms.awayArtifact,opts);

                end
            
            end
        end
    end

end

function [waveformsFiltered] = filterArtifactData(waveforms,opts)
    
    % pad
    waveformsFiltered = [waveforms,mean(waveforms(:,end-opts.PAD_MEAN_LENGTH:end),2)*ones(1,opts.NUM_PAD)]; %+ 1/2*randn(size(waveforms,1),opts.NUM_PAD).*sqrt(var(waveforms(:,end-opts.PAD_VAR_LENGTH:end),0,2))];
    % flip, filter and flip
    waveformsFiltered = fliplr(filter(opts.B,opts.A,fliplr(waveformsFiltered)')');
    
%     [b,a] = butter(2,7500/(30000/2),'low');
%     waveformsFiltered = fliplr(filter(b,a,fliplr(waveformsFiltered)')');
    % unpad
    waveformsFiltered = waveformsFiltered(:,1:end-opts.NUM_PAD);

end

function [figureHandle] = plotArtifacts(waveforms,opts)
    % deal with bit error
    if(opts.ADJUST_FOR_BIT_ERROR)
        waveforms = waveforms/0.254;
    end
    
    % deal with colors
    if(~isempty(opts.COLOR_RANGE))
        colorOrder = [linspace(opts.COLOR_RANGE(1,1),opts.COLOR_RANGE(2,1),size(waveforms,1))',...
            linspace(opts.COLOR_RANGE(1,2),opts.COLOR_RANGE(2,2),size(waveforms,1))',...
            linspace(opts.COLOR_RANGE(1,3),opts.COLOR_RANGE(2,3),size(waveforms,1))'];
    end

    % make figure
    figureHandle = figure();
%     xData  = ((1:size(waveforms,2))-1 - size(waveforms,2)/2)/30;
    xData  = ((1:size(waveforms,2))-1)/30;

    % enumerate over subplots
    for row = 1:opts.ROW_SUBPLOT
        for col = 1:opts.COL_SUBPLOT
            ax = subplot(opts.ROW_SUBPLOT,opts.COL_SUBPLOT,(row-1)*opts.COL_SUBPLOT + col);
            if(~isempty(opts.COLOR_RANGE))
                set(ax,'colorOrder',colorOrder);
                hold all
            end
            
            % plot
            plotIdx = ((row-1)*opts.COL_SUBPLOT + col-1)*opts.MAX_WAVES_PLOT + 1;
            if(plotIdx < size(waveforms,1))
                plot(xData,waveforms(plotIdx:min(size(waveforms,1),plotIdx+opts.MAX_WAVES_PLOT-1),:),'linewidth',1);
            end
            xlim(opts.XLIM)
            ylim(opts.YLIM)
            formatForLee(gcf)
        end
    end

   % title
    if(opts.PLOT_TITLE)
        suptitle(opts.TITLE_TO_PLOT)
    end
    
    % deal with saving plot
    if(opts.FIGURE_SAVE && strcmpi(opts.FIGURE_NAME,'')~=1 && strcmpi(opts.FIGURE_DIR,'')~=1)
        saveFiguresLIB(figHandle,opts.FIGURE_DIR,strcat(opts.FIGURE_NAME));
    end

end


function [opts] = configureOpts(optsInput)
    
    opts.ARTIFACT_MULTIPLIER = 1;
    opts.ADJUST_FOR_BIT_ERROR = 1;

    opts.MAX_WAVES_PLOT = 4;
    opts.ROW_SUBPLOT = 2;
    opts.COL_SUBPLOT = 2;
    
    opts.RANDOM_SAMPLE = 1;
    
    opts.TIME_AFTER_STIMULATION_ARTIFACT = 10/1000;
    
    opts.PLOT_FILTERED = [0,1];
    
    opts.PLOT_TITLE = 0; % not implemented
    opts.TITLE_TO_PLOT = '';
    
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = '';
    opts.FIGURE_PREFIX = '';
    
    opts.YLIM = [-500,500];
    opts.XLIM = [0,10];
    
    opts.COLOR_RANGE = []; %[0.0,0.0,0.0;0.5,0.5,0.5];
    
    [opts.B,opts.A] = butter(6,500/(30000/2),'high');
    opts.NUM_PAD = 200;
    opts.PAD_MEAN_LENGTH = 0;
    opts.PAD_VAR_LENGTH = 20;
    opts.WAVEFORMS_TO_PLOT = [];
    opts.CHANS_TO_PLOT = [];
    
    
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