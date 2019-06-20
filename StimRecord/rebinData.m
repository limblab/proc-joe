%% rebin

    optsExtract.BIN_SIZE = 1/1000; % s
    bin_edges = 1000*(-optsExtract.PRE_TIME:optsExtract.BIN_SIZE:optsExtract.POST_TIME);
    for u = 1:numel(arrayData)
        for cond = 1:numel(arrayData{u}.bC)
            arrayData{u}.bE{cond} = bin_edges;
            arrayData{u}.bC{cond} = histcounts(arrayData{u}.spikeTrialTimes{cond}*1000,bin_edges)/arrayData{u}.numStims(cond);
            arrayData{u}.binMaxYLim = max(arrayData{u}.binMaxYLim,1.1*max(arrayData{u}.bC{cond}));
        end
        
    end