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
    
    % get difference between pd from an electrode and the test angles
    pdHeatmapDiffData = zeros(size(pdHeatmapData,1),size(pdHeatmapData,2));
    pdHeatmapWeight = zeros(size(pdHeatmapData,1),size(pdHeatmapData,2));
    
    for angIdx = 1:numel(testAngles)
        for rowIdx = 1:size(pdHeatmapData,1)
            for colIdx = 1:size(pdHeatmapData,2)
                pdHeatmapDiffData(rowIdx,colIdx) = angleDiff(pdHeatmapData(rowIdx,colIdx),testAngles(angIdx),0,0);
%                 if(~isnan(stimHeatmapData.dataRatioScaled(rowIdx,colIdx)))
%                     weight_idx = find(pdHeatmapData(rowIdx,colIdx) > PDBinEdges,1,'last');
%                     pdHeatmapWeight(rowIdx,colIdx) = 1/PDBinCount(weight_idx);
%                 end
            end
        end
        % scale pd diff data to -1 -> 1
        pdHeatmapDiffData = -2*pdHeatmapDiffData/180 + 1;

        % multiply matrices
        matrixMultResult = stimHeatmapData.dataRatioScaled.*pdHeatmapDiffData;
        isEntry = ~isnan(stimHeatmapData.dataRatioScaled) & ~isnan(pdHeatmapDiffData);
        outputData.metricAllAngles(angIdx) = mean(matrixMultResult(isEntry));
    end
    
    [~,maxAngIdx] = max(outputData.metricAllAngles);
    outputData.bestAngle = testAngles(maxAngIdx);
    

end