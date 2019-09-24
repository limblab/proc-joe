function [arrayDataRebin] = rebinArrayData(arrayData,inputData)

    binSize = inputData.bin_size;

    
    
    arrayDataRebin = arrayData;
    bin_edges = (arrayData{1}.binEdges{1}(1)):binSize:(arrayData{1}.binEdges{1}(end)); % in ms
    for u = 1:numel(arrayDataRebin)
        for cond = 1:numel(arrayDataRebin{u}.binCounts)
            spike_trial_times = arrayDataRebin{u}.spikeTrialTimes{cond}; % in s for some reason
            
            
            arrayDataRebin{u}.binEdges{cond} = bin_edges;
            arrayDataRebin{u}.binCounts{cond} = histcounts(spike_trial_times*1000,bin_edges);
            
            if(max(arrayDataRebin{u}.binCounts{cond}) > arrayDataRebin{u}.binMaxYLim)
                arrayDataRebin{u}.binMaxYLim = max(arrayDataRebin{u}.binCounts{cond});
            end
            
        end
    end



end