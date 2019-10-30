function [figureHandles,heatmap_data] = plotHeatmaps(arrayData,mapFileName,opts)

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
    %% get heatmap data
    heatmap_data = {};
    if(opts.ALL_NEURONS)
        [heatmap_data] = getHeatmapDataAllNeurons(arrayData,opts);
    else
        heatmap_data = getHeatmapDataAllStimChans(arrayData,MAP_DATA,opts);
    end
    
    %% plot heatmaps
    for heatmap_idx = 1:numel(heatmap_data)
        dataRatio = heatmap_data{heatmap_idx}.dataRatioPlot;
        % plot heatmap
        figureHandles{end+1} = figure();
        figureHandles{end}.Position(4) = figureHandles{end}.Position(3);
        figureHandles{end}.Position(2) = figureHandles{end}.Position(2) - 200; % move down to not binEdges annoyingly off my screen

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
        % plot each data point
        for unit = 1:numel(dataRatio)
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

            rectangle('Position',[heatmap_data{heatmap_idx}.col(unit),(11-heatmap_data{heatmap_idx}.row(unit)),1,1],'EdgeColor','k',...
                        'FaceColor',colorToPlot,'linewidth',0.1);
                    
            plottedHere(heatmap_data{heatmap_idx}.col(unit),heatmap_data{heatmap_idx}.row(unit)) = 1;
%             unit_location(unit,:) = [heatmap_data{heatmap_idx}.col(unit),heatmap_data{heatmap_idx}.row(unit)];
%             unit_channel(unit) = heatmap_data{heatmap_idx}.chan(unit);
        end

        % plot stim chan
        posIdx = find(MAP_DATA.chan == heatmap_data{heatmap_idx}.main_chan);
        posList = [MAP_DATA.row,MAP_DATA.col];
        stimRow = posList(posIdx,1);
        stimCol = posList(posIdx,2);

        % magenta box
        plot([stimCol,stimCol+1,stimCol+1,stimCol,stimCol],[stimRow,stimRow,stimRow+1,stimRow+1,stimRow],'m','linewidth',3)

        if(~plottedHere(stimCol,stimRow))
            plot([stimCol,stimCol+1],[stimRow,stimRow+1],'m','linewidth',3);
            plot([stimCol+1,stimCol],[stimRow,stimRow+1],'m','linewidth',3);
        end

        set(gca,'Visible','off')
        axis square
%         unit_data{end+1}.data = dataRatio;
%         unit_data{end}.loc = unit_location;
%         unit_data{end}.chan = unit_channel;
        %% save figures
        if(opts.FIGURE_SAVE && strcmpi(opts.FIGURE_DIR,'')~=1)
            if(~opts.AUTO_WINDOW)
                FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_stimChan',num2str(heatmap_data{heatmap_idx}.main_chan),'_wave',num2str(heatmap_data{heatmap_idx}.wave),...
                    '_heatmap');
            else
                FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_stimChan',num2str(heatmap_data{heatmap_idx}.main_chan),'_wave',num2str(heatmap_data{heatmap_idx}.wave),...
                    '_heatmap');
            end
            saveFiguresLIB(figureHandles{end},opts.FIGURE_DIR,FIGURE_NAME);
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
                FIGURE_NAME = strcat(opts.FIGURE_PREFIX,'_heatmapBarPlot_',num2str(barMake));
                saveFiguresLIB(figureHandles{end},opts.FIGURE_DIR,FIGURE_NAME);
            end
        end
    end
    
end


function [heatmap_data] = getHeatmapDataAllNeurons(arrayData,opts)
    
    heatmap_idx = 1;
    for chan = opts.STIM_ELECTRODE_PLOT
        for wave = opts.WAVEFORM_TYPES_PLOT
           
            %define dataRatio
            dataRatio = [];

            for unit = 1:numel(arrayData)
                
                singleCondData = getHeatmapDataSingleCond(arrayData,unit,chan,wave,opts); 
                dataRatio(unit) = singleCondData.dataRatio;
                
                heatmap_data{heatmap_idx}.row(unit) = arrayData{unit}.ROW;
                heatmap_data{heatmap_idx}.col(unit) = arrayData{unit}.COL;
            end
            
            percentile = 90;
            dataRatioScaled = scaleDataRatio(dataRatio,percentile);

            if(opts.LOG_SCALE)
                dataRatio = dataRatio+eps;
                if(opts.MAX_RATIO == 0)
                    opts.MAX_RATIO = eps;
                end
            end

            if(opts.RELATIVE_INHIBITION)
                dataRatio(dataRatio < 0) = dataRatio(dataRatio < 0)./dataPre(dataRatio<0);
            end

            dataRatioPlot = dataRatio;
            dataRatioPlot(dataRatioPlot > opts.MAX_RATIO) = opts.MAX_RATIO;
            dataRatioPlot(dataRatioPlot < opts.MIN_RATIO) = opts.MIN_RATIO;
            
            heatmap_data{heatmap_idx}.dataRatio = dataRatio;
            heatmap_data{heatmap_idx}.dataRatioPlot = dataRatioPlot;
            heatmap_data{heatmap_idx}.dataRatioScaled = dataRatioScaled;
            heatmap_data{heatmap_idx}.chan = chan;
            heatmap_data{heatmap_idx}.wave = wave;
            heatmap_data{heatmap_idx}.main_chan = arrayData{1}.CHAN_LIST(chan);
            if(iscell(heatmap_data{heatmap_idx}.main_chan))
                heatmap_data{heatmap_idx}.main_chan = heatmap_data{heatmap_idx}.main_chan{1};
            end
            heatmap_idx = heatmap_idx + 1;
        end
    end

end


function [heatmap_data] = getHeatmapDataAllStimChans(arrayData,map_data,opts)

    heatmap_data = {};
    heatmap_idx = 1;
    for arrIdx = 1:numel(arrayData)
        for wave = opts.WAVEFORM_TYPES_PLOT
            dataPre = zeros(size(arrayData{arrIdx}.binCounts,1),1);
            dataPost = zeros(size(arrayData{arrIdx}.binCounts,1),1);

            for chan = opts.STIM_ELECTRODE_PLOT

                singleCondData = getHeatmapDataSingleCond(arrayData,arrIdx,chan,wave,opts); 
                dataRatio(chan) = singleCondData.dataRatio;

                map_data_idx = find(map_data.chan == arrayData{arrIdx}.CHAN_LIST{chan});

                heatmap_data{arrIdx}.row(chan) = 11-map_data.row(map_data_idx);
                heatmap_data{arrIdx}.col(chan) = map_data.col(map_data_idx);
            end

            
            percentile = 90;
            dataRatioScaled = scaleDataRatio(dataRatio,percentile);
            
            if(opts.LOG_SCALE)
                dataRatio = dataRatio+eps;
                if(opts.MAX_RATIO == 0)
                    opts.MAX_RATIO = eps;
                end
            end

            if(opts.RELATIVE_INHIBITION)
                dataRatio(dataRatio < 0) = dataRatio(dataRatio < 0)./dataPre(dataRatio<0);
            end

            dataRatioPlot = dataRatio;
            dataRatioPlot(dataRatioPlot > opts.MAX_RATIO) = opts.MAX_RATIO;
            dataRatioPlot(dataRatioPlot < opts.MIN_RATIO) = opts.MIN_RATIO;

            heatmap_data{heatmap_idx}.dataRatio = dataRatio;
            heatmap_data{heatmap_idx}.dataRatioPlot = dataRatioPlot;
            heatmap_data{heatmap_idx}.dataRatioScaled = dataRatioScaled;
            heatmap_data{heatmap_idx}.chan = 1:numel(arrayData{arrIdx}.CHAN_LIST);
            heatmap_data{heatmap_idx}.wave = wave*ones(numel(arrayData{arrIdx}.CHAN_LIST),1);
            heatmap_data{heatmap_idx}.main_chan = arrayData{arrIdx}.CHAN_REC;
            heatmap_idx = heatmap_idx + 1;
        end
    end

end

function [outputData] = getHeatmapDataSingleCond(arrayData,unit,chan,wave,opts)
    %initialize output data
    outputData = [];
    outputData.dataPre = 0; 
    outputData.dataPost = 0;
    outputData.numPostStimBins = 0;
    
    % get data in the correct ranges from arrayData
    postStim_binEdgePre = max(find(arrayData{unit}.binEdges{1} <= opts.STIM_PRE_TIME*1000));
    postStim_binEdgePost = max(find(arrayData{unit}.binEdges{1} <= opts.STIM_POST_TIME*1000));
    numPostStimBins = ones(numel(arrayData{unit}.binEdges),1)*(postStim_binEdgePost - postStim_binEdgePre);

    baseline_binEdgePre = max(find(arrayData{unit}.binEdges{1} <= opts.BASELINE_PRE_TIME*1000));
    baseline_binEdgePost = max(find(arrayData{unit}.binEdges{1} <= opts.BASELINE_POST_TIME*1000));
    
    try
        if(opts.AUTO_WINDOW && opts.EXCITATORY && arrayData{arrayDataIdx}.isExcitatory{chan,wave})
            tempPre = max(find(arrayData{unit}.binEdges{chan,wave} <= arrayData{unit}.excitatoryLatency{chan,wave}(1)));
            tempPost = max(find(arrayData{unit}.binEdges{chan,wave} <= arrayData{unit}.excitatoryLatency{chan,wave}(3)));

            dataPost = sum(arrayData{unit}.binCounts{chan,wave}(tempPre:tempPost));
            numPostStimBins = tempPost-tempPre;
        elseif(opts.AUTO_WINDOW && opts.INHIBITORY && arrayData{unit}.isInhibitory{chan,wave})
            tempPre = max(find(arrayData{unit}.binEdges{chan,wave} <= arrayData{unit}.inhibitoryLatency{chan,wave}(1)));
            tempPost = max(find(arrayData{unit}.binEdges{chan,wave} <= arrayData{unit}.inhibitoryLatency{chan,wave}(2)));

            dataPost = sum(arrayData{unit}.binCounts{chan,wave}(tempPre:tempPost));
            numPostStimBins = tempPost-tempPre;
        else
            dataPost = sum(arrayData{unit}.binCounts{chan,wave}(postStim_binEdgePre:postStim_binEdgePost));
        end
        dataPre = numPostStimBins(unit)*mean(arrayData{unit}.binCounts{chan,wave}(baseline_binEdgePre:baseline_binEdgePost));
    catch
        dataPre = 0;
        dataPost = 0;
        numPostStimBins = 0;
    end
    
     %defining arrays and bin sizess
    lastRep = arrayData{unit}.stimData{chan,wave}(numel(arrayData{unit}.stimData{chan,wave}));
    repList = cell(1,lastRep); %all spike times organized by repetition
    preRepList = zeros(1,lastRep); %number of baseline spikes per repetition
    postRepList = zeros(1,lastRep); %number of post-stim spikes per repetition

    %separating spike times by repetition into repList
    for spikeTime=1:numel(arrayData{unit}.spikeTrialTimes{chan,wave})
        rep = arrayData{unit}.stimData{chan,wave}(spikeTime);
        repList{1,rep} = [repList{1,rep}, arrayData{unit}.spikeTrialTimes{chan,wave}(spikeTime)];
    end

    %going through repList 
    %and counting baseline and post-stim spikes
    for rep=1:numel(repList)
        for spikeTime=1:numel(repList{1,rep})
            if (repList{1,rep}(spikeTime)) < opts.BASELINE_POST_TIME && (repList{1,rep}(spikeTime)) > opts.BASELINE_PRE_TIME 
                preRepList(rep) = preRepList(rep)+1;
            elseif (repList{1,rep}(spikeTime)) > opts.STIM_PRE_TIME && (repList{1,rep}(spikeTime)) < opts.STIM_POST_TIME
                postRepList(rep) = postRepList(rep)+1;
            end
        end
    end


    %calculating standardized value
    dataPreMean = mean(preRepList);
    dataPostMean = mean(postRepList);
    dataPreStd = std(preRepList);
    spikesStandardized = (dataPostMean-dataPreMean)/dataPreStd;
    
    % package outputs
    outputData.dataPreMean = dataPreMean; 
    outputData.dataPostMean = dataPostMean;
    outputData.dataPreStd = dataPreStd;
    outputData.dataRatio = spikesStandardized;
    outputData.numPostStimBins = numPostStimBins;
    
    
end

function [dataRatioScaled] = scaleDataRatio(dataRatio,percentile)

    % scale dataRatio to -1 to 1 for further analysis
    %find value of 90% percentile
    dataRatioScaled = dataRatio;
    maxValue = prctile(dataRatioScaled,percentile);

    for i=1:numel(dataRatioScaled)
        %if data point is above 90th percentile, set it to value of
        %90th percentile
        if dataRatioScaled(i)>maxValue
            dataRatioScaled(i)=maxValue;
        end

    end

    %scale from -1 to 1
    dataRatioScaled = 2*(dataRatioScaled-(min(dataRatioScaled)))/(maxValue-min(dataRatioScaled)) - 1;        


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
    
    opts.ALL_NEURONS = 1;
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