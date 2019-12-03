function [arrayDataRebin] = rebinArrayData(arrayData,binSize)
    
    remove_cell = 0;
    if(~iscell(arrayData))
        arrayData = {arrayData};
        remove_cell = 1;
    end
    % rebins data in arrayData based on inputData.bin_size
  
    
    arrayDataRebin = arrayData;
    binEdges = (arrayData{1}.binEdges{1}(1)):binSize:(arrayData{1}.binEdges{1}(end)); % in ms
    for u = 1:numel(arrayDataRebin)
        for cond = 1:numel(arrayDataRebin{u}.binCounts)
            spikeTrialTimes = arrayDataRebin{u}.spikeTrialTimes{cond}; % in s for some reason
            
            
            arrayDataRebin{u}.binEdges{cond} = binEdges;
            arrayDataRebin{u}.binCounts{cond} = histcounts(spikeTrialTimes*1000,binEdges);
            
            if(~isfield(arrayData{u},'binMaxYLim') || max(arrayDataRebin{u}.binCounts{cond}) > arrayDataRebin{u}.binMaxYLim)
                arrayDataRebin{u}.binMaxYLim = max(arrayDataRebin{u}.binCounts{cond});
            end
            
        end
    end

    if(remove_cell)
        arrayDataRebin = arrayDataRebin{1};
    end

end