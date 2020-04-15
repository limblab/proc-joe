function [outputData] = compareHeatmaps(stimHeatmapData,pdHeatmapData,testAngles,mapFileName)

    outputData.metricAllAngles = zeros(numel(testAngles),1);
    outputData.centerChanPD = [];
    outputData.bestAngle = [];
    
    mapData = loadMapFile(mapFileName);
    
    % compute pdHeatmapDiffData using stimHeatmapData.main_chan as the
    % relative channel
    
    mapDataIdx = find(mapData.chan == stimHeatmapData.main_chan);
    centerCol = mapData.col(mapDataIdx);
    centerRow = 11-mapData.row(mapDataIdx);
    stimHeatmapData.dataRatioScaled(centerRow,centerCol) = nan; % remove stimulated channel
    
    % output center chan PD
    outputData.centerChanPD = pdHeatmapData(centerRow,centerCol);
    
    pdHeatmapDiffData = zeros(size(pdHeatmapData,1),size(pdHeatmapData,2));
    
    for angIdx = 1:numel(testAngles)
        for rowIdx = 1:size(pdHeatmapData,1)
            for colIdx = 1:size(pdHeatmapData,2)
                pdHeatmapDiffData(rowIdx,colIdx) = angleDiff(pdHeatmapData(rowIdx,colIdx),testAngles(angIdx),0,0);
            end
        end
        % scale pd diff data to -1 -> 1
        pdHeatmapDiffData = -2*pdHeatmapDiffData/180 + 1;

        % multiple matrices
        matrixMultResult = stimHeatmapData.dataRatioScaled.*pdHeatmapDiffData;
        isEntry = ~isnan(stimHeatmapData.dataRatioScaled) & ~isnan(pdHeatmapDiffData);
        outputData.metricAllAngles(angIdx) = mean(matrixMultResult(isEntry));
    end
    
    [~,maxAngIdx] = 
    outputData.bestAngle = testAngles(maxAngIdx);
    

end