function [arrayDataRebin] = rebinArrayData(arrayData,binSize)
    
    remove_cell = 0;
    if(~iscell(arrayData))
        arrayData = {arrayData};
        remove_cell = 1;
    end
    % rebins data in arrayData based on inputData.bin_size
    arrayDataRebin = arrayData;
    
    
    
    for u = 1:numel(arrayDataRebin)
        % makes binEdges if it doesn't already exist
        if(~isfield(arrayData{u},'binEdges'))
            binLimits = round([min(cellfun(@min,arrayData{u}.spikeTrialTimes)),max(cellfun(@max,arrayData{u}.spikeTrialTimes))]*1000,-1); % round to nearest 10ms
            binEdges = binLimits(1):binSize:binLimits(2);
        else
            binEdges = (arrayData{u}.binEdges{1}(1)):binSize:(arrayData{u}.binEdges{1}(end)); % in ms
        end
        
        for condition = 1:numel(arrayDataRebin{u}.spikeTrialTimes)
            spikeTrialTimes = arrayDataRebin{u}.spikeTrialTimes{condition}; % in s for some reason
            
            
            arrayDataRebin{u}.binEdges{condition} = binEdges;
            arrayDataRebin{u}.binCounts{condition} = histcounts(spikeTrialTimes*1000,binEdges)./arrayDataRebin{u}.numStims(condition);
            
            if(~isfield(arrayData{u},'binMaxYLim') || max(arrayDataRebin{u}.binCounts{condition}) > arrayDataRebin{u}.binMaxYLim)
                arrayDataRebin{u}.binMaxYLim = max(arrayDataRebin{u}.binCounts{condition});
            end
            
        end
    end

    if(remove_cell)
        arrayDataRebin = arrayDataRebin{1};
    end

end