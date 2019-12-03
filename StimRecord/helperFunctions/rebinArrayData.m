function [arrayDataRebin] = rebinArrayData(arrayData,binSize)
    % rebins data in arrayData based on inputData.bin_size
  
    arrayDataRebin = arrayData;
    binEdges = (arrayData{1}.binEdges{1}(1)):binSize:(arrayData{1}.binEdges{1}(end)); % in ms
    for u = 1:numel(arrayDataRebin)
        for condition = 1:numel(arrayDataRebin{u}.binCounts)
            spikeTrialTimes = arrayDataRebin{u}.spikeTrialTimes{condition}; % in s for some reason
            
            
            arrayDataRebin{u}.binEdges{condition} = binEdges;
            arrayDataRebin{u}.binCounts{condition} = histcounts(spikeTrialTimes*1000,binEdges)./arrayDataRebin{u}.numStims(condition);
            
            if(~isfield(arrayData{u},'binMaxYLim') || max(arrayDataRebin{u}.binCounts{condition}) > arrayDataRebin{u}.binMaxYLim)
                arrayDataRebin{u}.binMaxYLim = max(arrayDataRebin{u}.binCounts{condition});
            end
            
        end
    end



end