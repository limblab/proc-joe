function [figureHandles] = plotHeatmaps(arrayData,mapFileName,opts)

    %% configure opts and set default values
    opts = configureOpts(opts);
    
    if(size(opts.STIM_ELECTRODE_PLOT,1) > size(opts.STIM_ELECTRODE_PLOT,2))
        opts.STIM_ELECTRODE_PLOT = opts.STIM_ELECTRODE_PLOT';
    end
    
    if(size(opts.WAVEFORM_TYPES_PLOT,1) > size(opts.WAVEFORM_TYPES_PLOT,2))
        opts.WAVEFORM_TYPES_PLOT = opts.WAVEFORM_TYPES_PLOT';
    end
    
    if(opts.LOG_SCALE)
        opts.MAX_RATIO = log10(opts.MAX_RATIO);
        opts.MIN_RATIO = log10(opts.MIN_RATIO);
    end
    
    %% useful constants
    figureHandles = {};
    
    NUM_CHANS = size(arrayData{1,1}.spikeTrialTimes,1);
    NUM_WAVEFORMS = size(arrayData{1,1}.spikeTrialTimes,2);
    
    MAP_DATA = loadMapFile(mapFileName);
    %% plot a heatmap for each waveform type and stimulation channel
    for chan = opts.STIM_ELECTRODE_PLOT
        for wave = opts.WAVEFORM_TYPES_PLOT
            %% get data in the correct ranges from arrayData
            postStim_binEdgePre = max(find(arrayData{1,1}.bE{chan,wave} <= opts.STIM_PRE_TIME*1000));
            postStim_binEdgePost = max(find(arrayData{1,1}.bE{chan,wave} <= opts.STIM_POST_TIME*1000)); - 1;
            
            baseline_binEdgePre = max(find(arrayData{1,1}.bE{chan,wave} <= opts.BASELINE_PRE_TIME*1000));
            baseline_binEdgePost = max(find(arrayData{1,1}.bE{chan,wave} <= opts.BASELINE_POST_TIME*1000)) - 1; % one more bin edge than bin count
            
            dataPre = zeros(size(arrayData,1),1);
            dataPost = zeros(size(arrayData,1),1);
            
            for unit = 1:size(arrayData,1)
                dataPost(unit) = mean(arrayData{unit,1}.bC{chan,wave}(postStim_binEdgePre:postStim_binEdgePost));
                dataPre(unit) = mean(arrayData{unit,1}.bC{chan,wave}(baseline_binEdgePre:baseline_binEdgePost));
            end
            
            dataRatio = dataPost./dataPre;
        
            if(opts.LOG_SCALE)
                dataRatio(dataRatio == 0) = 0+eps;
                dataRatio = log10(dataRatio);
            end
            
            dataRatio(dataRatio > opts.MAX_RATIO) = opts.MAX_RATIO;
            dataRatio(dataRatio < opts.MIN_RATIO) = opts.MIN_RATIO;
            
            
            %% plot heatmap
            figureHandles{end+1} = figure();
            figureHandles{end}.Position(4) = figureHandles{end}.Position(3);
            figureHandles{end}.Position(2) = figureHandles{end}.Position(2) - 200; % move down to not be annoyingly off my screen
            colors = colormap(jet);
            plottedHere = zeros(10,10);
            
            surface(zeros(opts.NUM_ROWS+1,opts.NUM_COLS+1));
            
            ax = gca;
            ax.Children(1).LineWidth = 1.5; % thicker box boundaries
            
            colormap([1,1,1]);
            
            ax.YTickLabel = {};
            ax.XTickLabel = {};
            
            hold on
            % plot each unit
            for unit = 1:size(arrayData,1)
                colorIdx = min(size(colors,1),max(1,floor(size(colors,1)*(dataRatio(unit)-opts.MIN_RATIO)/(opts.MAX_RATIO-opts.MIN_RATIO))));
                
                rectangle('Position',[arrayData{unit}.COL,arrayData{unit}.ROW,1,1],'EdgeColor','k',...
                            'FaceColor',colors(colorIdx,:),'linewidth',0.1);
                plottedHere(arrayData{unit}.COL,arrayData{unit}.ROW) = 1;
            end
            
            % plot stim chan
            posIdx = find(MAP_DATA.chan == arrayData{1,1}.CHAN_LIST(chan));
            posList = [MAP_DATA.row,MAP_DATA.col];
            stimRow = 11 - posList(posIdx,1);
            stimCol = posList(posIdx,2);
            
            % magenta box
            plot([stimCol,stimCol+1,stimCol+1,stimCol,stimCol],[stimRow,stimRow,stimRow+1,stimRow+1,stimRow],'m','linewidth',3)
            
            if(~plottedHere(stimCol,stimRow))
                plot([stimCol,stimCol+1],[stimRow,stimRow+1],'m','linewidth',3);
                plot([stimCol+1,stimCol],[stimRow,stimRow+1],'m','linewidth',3);
            end
            
            set(gca,'Visible','off')
            
            
            %% save figures
            if(opts.FIGURE_SAVE && strcmpi(FIGURE_DIR,'')~=1)
                FIGURE_NAME = strcat(FIGURE_PREFIX,'_stimChan',num2str(arrayData{1,1}.CHAN_LIST(chan)),'_wave',num2str(wave),'_heatmap');
                saveFiguresLIB(figHandle,optsSave.FIGURE_DIR,FIGURE_NAME);
            end
            
            
        end
    end
    
    %% make bar plot
    if(opts.MAKE_BAR_PLOT)
        figureHandles{end+1} = figure();
        b = colorbar;
        colormap jet;
        set(gca,'Visible','off');
        b.FontSize = 14;
        b.Ticks = [0,0.25,0.5,0.75,1.0];
        
        maxDataRound = round(opts.MAX_RATIO,1);
        minDataRound = round(opts.MIN_RATIO,1);
        
        b.TickLabels = {};
        for i = b.Ticks
            if(i==b.Ticks(end))
                b.TickLabels{end+1,1} = strcat('>',num2str(i*(maxDataRound-minDataRound) + minDataRound));
            else
                b.TickLabels{end+1,1} = num2str(i*(maxDataRound-minDataRound) + minDataRound);
            end
        end
        
        if(opts.FIGURE_SAVE && strcmpi(FIGURE_DIR,'')~=1)
            FIGURE_NAME = strcat(FIGURE_PREFIX,'_stimChan',num2str(arrayData{1,1}.CHAN_LIST(chan)),'_wave',num2str(wave),'_heatmapBarPlot');
            saveFiguresLIB(figHandle,optsSave.FIGURE_DIR,FIGURE_NAME);
        end
    end
    
end



function [opts] = configureOpts(optsInput)

    opts.MAKE_BAR_PLOT = 1;

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
    
    opts.NUM_ROWS = 10;
    opts.NUM_COLS = 10;
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