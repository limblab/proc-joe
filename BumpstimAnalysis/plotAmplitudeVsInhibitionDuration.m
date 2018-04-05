function [figureHandles,FITS] = plotAmplitudeVsInhibitionDuration(arrayData,mapFileName,opts)
    
    %% configure opts and set default values
    opts = configureOpts(opts);
    
    if(size(opts.STIM_ELECTRODE_PLOT,1) > size(opts.STIM_ELECTRODE_PLOT,2))
        opts.STIM_ELECTRODE_PLOT = opts.STIM_ELECTRODE_PLOT';
    end
    
    if(size(opts.WAVEFORM_TYPES_PLOT,1) > size(opts.WAVEFORM_TYPES_PLOT,2))
        opts.WAVEFORM_TYPES_PLOT = opts.WAVEFORM_TYPES_PLOT';
    end
    
    %% useful constants
    figureHandles = {};
    FITS = [];
    
    NUM_CHANS = size(arrayData{1,1}.spikeTrialTimes,1);
    NUM_WAVEFORMS = size(arrayData{1,1}.spikeTrialTimes,2);
    
    MAP_DATA = loadMapFile(mapFileName);

    
    if(opts.PLOT_ON_ONE_FIGURE)
        figureHandles{end+1} = figure();
        hold on
    end

    for chan = opts.STIM_ELECTRODE_PLOT
        for wave = opts.WAVEFORM_TYPES_PLOT
            %% get data in the correct ranges from arrayData
            postStim_binEdgePre = max(find(arrayData{1,1}.bE{chan,wave} <= opts.STIM_PRE_TIME*1000));
            postStim_binEdgePost = max(find(arrayData{1,1}.bE{chan,wave} <= opts.STIM_POST_TIME*1000)) - 1;
            numPostStimBins = postStim_binEdgePost - postStim_binEdgePre;
            
            baseline_binEdgePre = max(find(arrayData{1,1}.bE{chan,wave} <= opts.BASELINE_PRE_TIME*1000));
            baseline_binEdgePost = max(find(arrayData{1,1}.bE{chan,wave} <= opts.BASELINE_POST_TIME*1000)) - 1; % one more bin edge than bin count
            
            dataPre = zeros(size(arrayData,1),1);
            dataPost = zeros(size(arrayData,1),1);
            inhibitionDuration = zeros(size(arrayData,1),1);
            
            for unit = 1:size(arrayData,1)
                dataPost(unit) = sum(arrayData{unit,1}.bC{chan,wave}(postStim_binEdgePre:postStim_binEdgePost));
                dataPre(unit) = numPostStimBins*mean(arrayData{unit,1}.bC{chan,wave}(baseline_binEdgePre:baseline_binEdgePost));
                inhibitionDuration(unit) = diff(arrayData{unit}.inhibitoryLatency{chan,wave});
            end
            
            dataRatio = dataPost-dataPre;
            
            dataRatio(inhibitionDuration <= 0) = [];
            inhibitionDuration(inhibitionDuration <= 0) = [];
            
            %% plot distance vs amplitude
            if(~opts.PLOT_ON_ONE_FIGURE)
                figureHandles{end+1} = figure();
                c = 'k';
            else
                c = opts.COLORS{chan*wave};
            end
            plot(inhibitionDuration,dataRatio,'.','markersize',opts.MARKER_SIZE,'color',c);
%             hold on
%             [FITS.f{chan*wave},FITS.stats{chan*wave}] = fit(distances,dataRatio,'exp1');
%             plot(FITS.f{chan*wave});
            %% save figures
            if(opts.FIGURE_SAVE && strcmpi(FIGURE_DIR,'')~=1)
                FIGURE_NAME = strcat(FIGURE_PREFIX,'_stimChan',num2str(arrayData{1,1}.CHAN_LIST(chan)),'_wave',num2str(wave),'_heatmap');
                saveFiguresLIB(figHandle,optsSave.FIGURE_DIR,FIGURE_NAME);
            end
            
            
        end
    end
    
    
end



function [opts] = configureOpts(optsInput)

    opts.BASELINE_PRE_TIME = -10/1000;
    opts.BASELINE_POST_TIME = -3/1000;
    opts.STIM_PRE_TIME = 1.2/1000;
    opts.STIM_POST_TIME = 5/1000;
    
    opts.STIM_ELECTRODE_PLOT = 1;
    opts.WAVEFORM_TYPES_PLOT = 1;
    
    opts.MARKER_SIZE = 20;
    
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = '';
    opts.FIGURE_PREFIX = '';
    
    opts.PLOT_ON_ONE_FIGURE = 1;
    opts.COLORS = {'r','b','g','k','m',[0.5,0.5,0.5],'y'};

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