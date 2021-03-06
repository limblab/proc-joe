function [outputData] = combineHeatmaps(inputData)

%need dataRatioScaled (spikes) -- array (1x23) -- from analyzeStimData
%need PDscaled (preferred directions) -- array (10x10) -- from analyzeCObump
%also need... td_all -- from analyzeCObump, 
%             mapData -- from analyzeCObump, 
%             arrayData -- from analyzeStimData.

    % make some output variables
    figure_handle = [];
    heatmap = [];
    
    %extract variables from input_data
    mapData = inputData.mapData;
    arrayData = inputData.arrayData;
    dataRatioScaled = inputData.dataRatioScaled;
    PDscaled = inputData.PDscaled;
    waveform = inputData.waveform;
    angleNum = inputData.angleNumber;
    mainChannel = inputData.mainChan;

    %other variables needed
    plotChannels = zeros(1,numel(arrayData));
    for i=1:numel(arrayData)
        plotChannels(i) = arrayData{1,i}.CHAN_REC; %list of channels with units on them
    end
    mainRow = 11-mapData.row(find(mapData.chan == mainChannel)); %stimchannel row
    mainCol = mapData.col(find(mapData.chan == mainChannel)); %stimchannel column

    if inputData.waveform == 0
        figureName = strcat(inputData.figPrefix,'_mainChan',num2str(mainChannel),'_angleNumber',num2str(angleNum),'_combinedHeatmap');
    else
        figureName = strcat(inputData.figPrefix,'_mainChan',num2str(mainChannel),'_waveForm',num2str(waveform),'_combinedHeatmap');
    end
    alphaArray = zeros(10,10);
    
    %create 10x10 matrix version of dataRatioScaled
    dataRatioMatrix = nan(10,10);
    for i=1:numel(dataRatioScaled)
        index = find(mapData.chan == plotChannels(i));
        dataRatioMatrix(11-mapData.row(index),mapData.col(index)) = dataRatioScaled(i);
        alphaArray(11-mapData.row(index),mapData.col(index)) = 1;
    end

    %multiply arrays together
    heatmapData = PDscaled.*dataRatioMatrix;

%     disp(heatmapData)

    %plot heatmap
    if inputData.plotHeatmap == 1
        [figure_handle,heatmap] = plotHeatmap(heatmapData,alphaArray,figureName,mainRow,mainCol);
    end
    
    %calculate average
    average = sum(heatmapData(~isnan(heatmapData)))/sum(sum(~isnan(heatmapData)));
    
    % package outputs
    outputData.heatmapData = heatmapData;
    outputData.figure_handle = figure_handle;
    outputData.heatmap = heatmap;
    outputData.average = average;
    
end

%plot heatmap
function [f,heatmap] = plotHeatmap(heatmapData,alphaArray,figureName,mainRow,mainCol)
    %plot the pds
    f = figure();
    f.Name = figureName;
    clims=[-1,1];
    heatmap = imagesc(heatmapData,'alphaData',alphaArray,clims);

    colormap((viridis));
    colorbar;

    % magenta box for main chan
    hold on
    rectangle('position',[(mainCol-0.5) (mainRow-0.5) 1 1],'edgecolor','m','linewidth',2);
    plot([mainCol-0.5,mainCol+0.5],[mainRow-0.5,mainRow+0.5],'m','linewidth',2);
    plot([mainCol+0.5,mainCol-0.5],[mainRow-0.5,mainRow+0.5],'m','linewidth',2);

end

%plot average combinedData value vs. center angle (1-360)
function [] = plotCombinedMeanValues()
end