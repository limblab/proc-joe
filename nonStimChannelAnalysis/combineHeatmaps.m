function [outputData] = combineHeatmaps(inputData)

%need dataRatioScaled (spikes) -- array (1x23) -- from analyzeStimData
%need PDscaled (preferred directions) -- array (96x1) -- from analyzeCObump
%also need... td_all -- from analyzeCObump, 
%             mapData -- from analyzeCObump, 
%             arrayData -- from analyzeStimData.

    % extract variables from input_data
    td_all = inputData.td_all;
    mapData = inputData.mapData;
    arrayData = inputData.arrayData;
    dataRatioScaled = inputData.dataRatioScaled;
    PDscaled = input_data.PDscaled;

    %other variables needed
    plotChannels = zeros(1,numel(arrayData));
    for i=1:numel(arrayData)
        plotChannels(i) = arrayData{1,i}.CHAN_REC; %list of channels with units on them
    end
    stimChannel = 9; %channel being stimulated
    figPrefix = 'Han_20190924'; %prefix for figure name
    figureName = strcat(figPrefix,'_stimChan',num2str(stimChannel),'_combinedHeatmap');
    heatmapData = nan(10,10);
    alphaArray = zeros(10,10);

    %create new PDscaled array of just channels with units -- array (1x23)
    PDscaledUnits = [];
    for i=plotChannels
        PDscaledUnits = [PDscaledUnits PDscaled(i)];
    end

    %multiply arrays together
    combinedData = PDscaledUnits.*dataRatioScaled;

    %put channel values in correct location in new array
    for pd=plotChannels
        for chan=1:96
            if td_all(1).LeftS1_unit_guide(pd) == mapData.chan(chan)
                disp(pd)
                disp(chan)
                heatmapData((11-mapData.row(chan)),mapData.col(chan)) = combinedData(find(plotChannels==pd));
                alphaArray(11-mapData.row(chan),mapData.col(chan))=1;
                %setting stimrow and stimcol
                if mapData.chan(chan) == stimChannel
                    stimRow = 11 - mapData.row(chan);
                    stimCol = mapData.col(chan);
                end
            end
        end
    end

    disp(heatmapData)

    %plot heatmap
    [figure_handle,heatmap] = plotHeatmap(heatmapData,alphaArray,figureName,stimRow,stimCol);
    
    
    % package outputs
    outputData.heatmapData = heatmapData;
    outputData.figure_handle = figure_handle;
    outputData.heatmap = heatmap;
    
end


function [f,heatmap] = plotHeatmap(heatmapData,alphaArray,figureName,stimRow,stimCol)
    %plot the pds
    f = figure();
    f.Name = figureName;
    clims=[-1,1];
    heatmap = imagesc(heatmapData,'alphaData',alphaArray,clims);

    colormap(flip(viridis,1));
    colorbar;

    % magenta box for stim chan
    hold on
    rectangle('position',[(stimCol-0.5) (stimRow-0.5) 1 1],'edgecolor','m','linewidth',2);
    plot([stimCol-0.5,stimCol+0.5],[stimRow-0.5,stimRow+0.5],'m','linewidth',2);
    plot([stimCol+0.5,stimCol-0.5],[stimRow-0.5,stimRow+0.5],'m','linewidth',2);

end