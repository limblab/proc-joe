function [figureHandles,FITS] = plotAmplitudeVsDistance(arrayData,mapFileName,opts)
    
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
            postStim_binEdgePre = max(find(arrayData{1}.bE{chan,wave} <= opts.STIM_PRE_TIME*1000));
            postStim_binEdgePost = max(find(arrayData{1}.bE{chan,wave} <= opts.STIM_POST_TIME*1000)) - 1;
            numPostStimBins = ones(numel(arrayData),1)*(postStim_binEdgePost - postStim_binEdgePre);
            
            baseline_binEdgePre = max(find(arrayData{1}.bE{chan,wave} <= opts.BASELINE_PRE_TIME*1000));
            baseline_binEdgePost = max(find(arrayData{1}.bE{chan,wave} <= opts.BASELINE_POST_TIME*1000)) - 1; % one more bin edge than bin count
            
            dataPre = zeros(numel(arrayData),1);
            dataPost = zeros(numel(arrayData),1);
            distances = zeros(numel(arrayData),1);
            STIMCHAN_POS = [11-MAP_DATA.row(find(MAP_DATA.chan == arrayData{1}.CHAN_LIST(chan))), MAP_DATA.col(find(MAP_DATA.chan == arrayData{1}.CHAN_LIST(chan)))];
            
            for unit = 1:numel(arrayData)
                if(opts.AUTO_WINDOW && arrayData{unit}.isExcitatory{chan,wave})
                    tempPre = max(find(arrayData{1}.bE{chan,wave} <= arrayData{unit}.excitatoryLatency{chan,wave}(1)));
                    tempPost = max(find(arrayData{1}.bE{chan,wave} <= arrayData{unit}.excitatoryLatency{chan,wave}(3)));
                    
                    dataPost(unit) = sum(arrayData{unit}.bC{chan,wave}(tempPre:tempPost));
                    numPostStimBins(unit) = tempPost-tempPre;
                else
                    dataPost(unit) = sum(arrayData{unit}.bC{chan,wave}(postStim_binEdgePre:postStim_binEdgePost));
                end
                dataPre(unit) = numPostStimBins(unit)*mean(arrayData{unit}.bC{chan,wave}(baseline_binEdgePre:baseline_binEdgePost));

                distances(unit) = 400*sqrt((arrayData{unit}.ROW-STIMCHAN_POS(1)).^2 + (arrayData{unit}.COL-STIMCHAN_POS(2)).^2);
            end
            
            dataRatio = dataPost-dataPre;
            
            %% plot distance vs amplitude
            if(~opts.PLOT_ON_ONE_FIGURE)
                figureHandles{end+1} = figure();
                c = 'k';
            else
                c = opts.COLORS{chan*wave};
            end
            plot(distances,dataRatio,'.','markersize',opts.MARKER_SIZE,'color',c);
            ylim([0,1])
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
    opts.COLORS = {'k','r','b',[0,0.5,0],'m'};

    opts.AUTO_WINDOW = 0;
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