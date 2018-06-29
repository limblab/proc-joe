function [figureHandles,unit_data] = plotHeatmaps(arrayData,mapFileName,opts)

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
    
    NUM_CHANS = size(arrayData{1,1}.spikeTrialTimes,1);
    NUM_WAVEFORMS = size(arrayData{1,1}.spikeTrialTimes,2);
    
    MAP_DATA = loadMapFile(mapFileName);
    
    colorsBelowOne = [0*(0:63)' 0*(0:63)' 0+(0:63)'/(63)];
    colorsAboveOne = [0+(0:63)'/(63) 0*(0:63)' 0*(0:63)'];
    %% plot a heatmap for each waveform type and stimulation channel
    for chan = opts.STIM_ELECTRODE_PLOT
        for wave = opts.WAVEFORM_TYPES_PLOT
            %% get data in the correct ranges from arrayData
            postStim_binEdgePre = max(find(arrayData{1}.bE{chan,wave} <= opts.STIM_PRE_TIME*1000));
            postStim_binEdgePost = max(find(arrayData{1}.bE{chan,wave} <= opts.STIM_POST_TIME*1000));
            numPostStimBins = ones(numel(arrayData))*(postStim_binEdgePost - postStim_binEdgePre);

            baseline_binEdgePre = max(find(arrayData{1}.bE{chan,wave} <= opts.BASELINE_PRE_TIME*1000));
            baseline_binEdgePost = max(find(arrayData{1}.bE{chan,wave} <= opts.BASELINE_POST_TIME*1000));

            dataPre = zeros(numel(arrayData),1);
            dataPost = zeros(numel(arrayData),1);

            for unit = 1:numel(arrayData)
                try
                    if(opts.AUTO_WINDOW && opts.EXCITATORY && arrayData{unit}.isExcitatory{chan,wave})
                        tempPre = max(find(arrayData{1}.bE{chan,wave} <= arrayData{unit}.excitatoryLatency{chan,wave}(1)));
                        tempPost = max(find(arrayData{1}.bE{chan,wave} <= arrayData{unit}.excitatoryLatency{chan,wave}(3)));

                        dataPost(unit) = sum(arrayData{unit}.bC{chan,wave}(tempPre:tempPost));
                        numPostStimBins(unit) = tempPost-tempPre;
                    elseif(opts.AUTO_WINDOW && opts.INHIBITORY && arrayData{unit}.isInhibitory{chan,wave})
                        tempPre = max(find(arrayData{1}.bE{chan,wave} <= arrayData{unit}.inhibitoryLatency{chan,wave}(1)));
                        tempPost = max(find(arrayData{1}.bE{chan,wave} <= arrayData{unit}.inhibitoryLatency{chan,wave}(2)));

                        dataPost(unit) = sum(arrayData{unit}.bC{chan,wave}(tempPre:tempPost));
                        numPostStimBins(unit) = tempPost-tempPre;
                    else
                        dataPost(unit) = sum(arrayData{unit}.bC{chan,wave}(postStim_binEdgePre:postStim_binEdgePost));
                    end
                    dataPre(unit) = numPostStimBins(unit)*mean(arrayData{unit}.bC{chan,wave}(baseline_binEdgePre:baseline_binEdgePost));
                catch
                    dataPre(unit) = 0;
                    dataPost(unit) = 0;
                    numPostStimBins(unit) = 0;
                end
            end

            dataRatio = dataPost-dataPre;

            if(opts.LOG_SCALE)
                dataRatio = dataRatio+eps;
                if(opts.MAX_RATIO == 0)
                    opts.MAX_RATIO = eps;
                end
            end

            if(opts.RELATIVE_INHIBITION)
                dataRatio(dataRatio < 0) = dataRatio(dataRatio < 0)./dataPre(dataRatio<0);
            end

            dataRatio(dataRatio > opts.MAX_RATIO) = opts.MAX_RATIO;
            dataRatio(dataRatio < opts.MIN_RATIO) = opts.MIN_RATIO;
            

            %% plot heatmap
            figureHandles{end+1} = figure();
            figureHandles{end}.Position(4) = figureHandles{end}.Position(3);
            figureHandles{end}.Position(2) = figureHandles{end}.Position(2) - 200; % move down to not be annoyingly off my screen

            plottedHere = zeros(10,10);

            % no grid lines
%             surface(zeros(opts.NUM_ROWS+1,opts.NUM_COLS+1));
            plot([1,opts.NUM_ROWS+1,opts.NUM_ROWS+1,1,1],[1,1,opts.NUM_COLS+1,opts.NUM_COLS+1,1],'-k','linewidth',1.5)

            ax = gca;
%             ax.Children(1).LineWidth = 1.5; % thicker box boundaries

            colormap([1,1,1]);

            ax.YTickLabel = {};
            ax.XTickLabel = {};

            hold on
            % plot each unit
            for unit = 1:numel(arrayData)
                if(dataRatio(unit) < 0)
                    if(~opts.LOG_SCALE)
                        mapping = floor(size(colorsBelowOne,1)*(-dataRatio(unit))/(-opts.MIN_RATIO));
                    else
                        mapping = floor(size(colorsBelowOne,1)*(log10(1+-(opts.LOG_PARAM-1)*dataRatio(unit)/-opts.MIN_RATIO)/log10(opts.LOG_PARAM)));
%                         mapping = floor(size(colorsBelowOne,1)*(1-(log10(-dataRatio(unit)/-opts.MIN_RATIO))/(log10(eps))));
                    end
                    colorIdx = min(size(colorsBelowOne,1),max(1,mapping));
                    colorToPlot = colorsBelowOne(colorIdx,:);
                elseif(dataRatio(unit) >= 0)
                    if(~opts.LOG_SCALE)
                        mapping = floor(size(colorsAboveOne,1)*(dataRatio(unit))/(opts.MAX_RATIO));
                    else
                        mapping = floor(size(colorsAboveOne,1)*(log10(1 + opts.LOG_PARAM*dataRatio(unit)/opts.MAX_RATIO)/log10(1+opts.LOG_PARAM)));
%                         mapping = floor(size(colorsBelowOne,1)*(1-(log10(dataRatio(unit)/opts.MAX_RATIO))/(log10(eps))));
                    end
                    colorIdx = min(size(colorsAboveOne,1),max(1,mapping));
                    colorToPlot = colorsAboveOne(colorIdx,:);
                end

                rectangle('Position',[arrayData{unit}.COL,arrayData{unit}.ROW,1,1],'EdgeColor','k',...
                            'FaceColor',colorToPlot,'linewidth',0.1);
                plottedHere(arrayData{unit}.COL,arrayData{unit}.ROW) = 1;
                unit_location(unit,:) = [arrayData{unit}.COL,arrayData{unit}.ROW];
                unit_channel(unit) = arrayData{unit}.CHAN_REC;
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
            axis square
            unit_data{chan,wave}.data = dataRatio;
            unit_data{chan,wave}.loc = unit_location;
            unit_data{chan,wave}.chan = unit_channel;
            %% save figures
            if(opts.FIGURE_SAVE && strcmpi(opts.FIGURE_DIR,'')~=1)
                if(~opts.AUTO_WINDOW)
                    FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_stimChan',num2str(arrayData{1,1}.CHAN_LIST(chan)),'_wave',num2str(wave),...
                        '_window',num2str(opts.STIM_PRE_TIME),'-',num2str(opts.STIM_POST_TIME),'_relativeInhib',num2str(opts.RELATIVE_INHIBITION),'_heatmap');
                else
                    FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_stimChan',num2str(arrayData{1,1}.CHAN_LIST(chan)),'_wave',num2str(wave),...
                        '_windowAUTO','_relativeInhib',num2str(opts.RELATIVE_INHIBITION),'_heatmap');
                end
                saveFiguresLIB(figureHandles{end},opts.FIGURE_DIR,FIGURE_NAME);
            end

        end
    end
    
    %% make bar plots
    if(opts.MAKE_BAR_PLOT)
        for barMake = 1:2
            figureHandles{end+1} = figure();
            b = colorbar;
            if(barMake == 1)
                colormap(colorsAboveOne);
            else
                colormap(flip(colorsBelowOne,1));
            end
            set(gca,'Visible','off');
            b.FontSize = 14;
            b.Ticks = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]; 
            b.TickDirection = 'out';
            
            if(opts.LOG_SCALE)
                b.Ticks = log10(1+(opts.LOG_PARAM-1)*b.Ticks)/log10(opts.LOG_PARAM);
%                 if(barMake == 1)
%                     b.Ticks = 1-(log10((b.Ticks+eps)*opts.MAX_RATIO))/(log10(eps));  
%                 else
%                     b.Ticks = 1-(log10((b.Ticks+eps)*opts.MIN_RATIO))/(log10(eps));
%                 end
            end
            
            if(barMake == 1)
                maxDataRound = round(opts.MAX_RATIO,1);
                minDataRound = 0;
            else
                maxDataRound = 0;
                minDataRound = round(opts.MIN_RATIO,1);
            end
            
            b.TickLabels = cell(1,numel(b.Ticks));

            for i = 1:2:numel(b.Ticks)
                if(i==numel(b.Ticks))
                    b.TickLabels{i,1} = strcat('>',num2str((i-1)*(maxDataRound-minDataRound)/(numel(b.Ticks)-1) + minDataRound));
                elseif(i==1)
                    b.TickLabels{i,1} = strcat('<',num2str((i-1)*(maxDataRound-minDataRound)/(numel(b.Ticks)-1) + minDataRound));
                else
                    b.TickLabels{i,1} = num2str((i-1)*(maxDataRound-minDataRound)/(numel(b.Ticks)-1) + minDataRound);
                end

            end

            if(opts.FIGURE_SAVE && strcmpi(opts.FIGURE_DIR,'')~=1)
                FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_stimChan',num2str(arrayData{1,1}.CHAN_LIST(chan)),'_wave',num2str(wave),...
                    '_window',num2str(opts.STIM_PRE_TIME),'-',num2str(opts.STIM_POST_TIME),'_relativeInhib',num2str(opts.RELATIVE_INHIBITION),...
                    '_heatmapBarPlot_',num2str(barMake));
                saveFiguresLIB(figureHandles{end},opts.FIGURE_DIR,FIGURE_NAME);
            end
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
    
    opts.MAX_RATIO = 1;
    opts.MIN_RATIO = -0.2;
    opts.LOG_SCALE = 0;
    opts.LOG_PARAM = 9;
    
    opts.FIGURE_SAVE = 0;
    opts.FIGURE_DIR = '';
    opts.FIGURE_PREFIX = '';
    
    opts.NUM_ROWS = 10;
    opts.NUM_COLS = 10;
    
    opts.RELATIVE_INHIBITION = 0;
    opts.AUTO_WINDOW = 0;
    opts.INHIBITORY = 0;
    opts.EXCITATORY = 1;
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