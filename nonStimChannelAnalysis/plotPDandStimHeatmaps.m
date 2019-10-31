function [outputData] = plotPDandStimHeatmaps(stimHeatmapData,pdHeatmapData,mapFileName)

    outputData = [];
    mapData = loadMapFile(mapFileName);
    
    % compute pdHeatmapDiffData using stimHeatmapData.main_chan as the
    % relative channel
    
    mapDataIdx = find(mapData.chan == stimHeatmapData.main_chan);
    centerCol = mapData.col(mapDataIdx);
    centerRow = 11-mapData.row(mapDataIdx);
    
    PDCenter = pdHeatmapData(centerRow,centerCol);
    
    for rowIdx = 1:size(pdHeatmapData,1)
        for colIdx = 1:size(pdHeatmapData,2)
            pdHeatmapDiffData(rowIdx,colIdx) = angleDiff(pdHeatmapData(rowIdx,colIdx),PDCenter,0,0);
        end
    end
    
    pdHeatmapScaled = -2*pdHeatmapDiffData/180 + 1;
    
    % alpha array is all spots that have entries in both matrices
    alphaArray = nan(size(stimHeatmapData.dataRatioScaled));
    alphaArray(~isnan(stimHeatmapData.dataRatioScaled) & ~isnan(pdHeatmapDiffData)) = 1;
    
    outputData.figure_handle = figure();
    outputData.figure_handle.Position = [202 440 1669 420];
    
    % plot PD difference heatmap
    subplot(1,3,1)
    PDHeatmap = imagesc(pdHeatmapDiffData,'alphaData',alphaArray);
    axis square
    colormap(gca,flip(viridis,1));
    caxis([0,180])
    b=colorbar; b.Limits = [0,180];
    % magenta box for stim chan
    hold on
    rectangle('position',[(centerCol-0.5) (centerRow-0.5) 1 1],'edgecolor','m','linewidth',2);
    plot([centerCol-0.5,centerCol+0.5],[centerRow-0.5,centerRow+0.5],'m','linewidth',2);
    plot([centerCol+0.5,centerCol-0.5],[centerRow-0.5,centerRow+0.5],'m','linewidth',2);
    
    % plot stim heatmap
    subplot(1,3,2)
    stimHeatmap = imagesc(stimHeatmapData.dataRatioScaled,'alphaData',alphaArray);
    axis square
    colormap(gca,viridis);
    caxis([-1,1])
    b=colorbar; b.Limits = [-1,1];
    % magenta box for stim chan
    hold on
    rectangle('position',[(centerCol-0.5) (centerRow-0.5) 1 1],'edgecolor','m','linewidth',2);
    plot([centerCol-0.5,centerCol+0.5],[centerRow-0.5,centerRow+0.5],'m','linewidth',2);
    plot([centerCol+0.5,centerCol-0.5],[centerRow-0.5,centerRow+0.5],'m','linewidth',2);
    
    % plot multiplied matrix
    subplot(1,3,3)
    stimHeatmap = imagesc(stimHeatmapData.dataRatioScaled.*pdHeatmapScaled,'alphaData',alphaArray);
    axis square
    colormap(gca,viridis);
    caxis([-1,1])
    b=colorbar; b.Limits = [-1,1];
    % magenta box for stim chan
    hold on
    rectangle('position',[(centerCol-0.5) (centerRow-0.5) 1 1],'edgecolor','m','linewidth',2);
    plot([centerCol-0.5,centerCol+0.5],[centerRow-0.5,centerRow+0.5],'m','linewidth',2);
    plot([centerCol+0.5,centerCol-0.5],[centerRow-0.5,centerRow+0.5],'m','linewidth',2);
    

end