function [arrayData] = getProbabilityOfResponse(arrayData,opts)

    opts = configureOpts(opts);
    binSize = mean(diff(arrayData{1}.bE{1,1})); % in ms
    %% calculate probability of detection for each channel
    for arrayDataIdx = 1:numel(arrayData)
        meanSpikesEvoked_all = zeros(size(arrayData{arrayDataIdx}.bC));
        normNumStimsResponsive = zeros(size(arrayData{arrayDataIdx}.bC));
        meanSpikesEvoked_responsive = zeros(size(arrayData{arrayDataIdx}.bC));
        
        
        for c = 1:size(arrayData{arrayDataIdx}.bE,1)
            for w = 1:1:size(arrayData{arrayDataIdx}.bE,2)
                % determine window if automatic windowing
                if(opts.AUTOMATIC_WINDOW) % get high time and low time automatically
                    opts.WINDOW = getAutomaticWindow(arryaData,arrayDataIdx,c,w,opts);
                end
                
                [meanSpikesEvoked_all(c,w), normNumStimsResponsive(c,w), meanSpikesEvoked_responsive(c,w)] = computeProbabilityOfResponse(arrayData,arrayDataIdx,c,w,opts);
                
            end
        end
        
        arrayData{arrayDataIdx}.singleProb.meanSpikesEvoked_all = meanSpikesEvoked_all;
        arrayData{arrayDataIdx}.singleProb.normNumStimsResponsive = normNumStimsResponsive;
        arrayData{arrayDataIdx}.singleProb.meanSpikesEvoked_responsive = meanSpikesEvoked_responsive;

    end
    
    %% compute pairwise probabilities
    
    for arrayDataIdx = 1:numel(arrayData)
        arrayData{arrayDataIdx}.pairwiseProb.normNumStimsWithBothResponse = zeros(numel(arrayData),size(arrayData{arrayDataIdx}.bC,1),size(arrayData{arrayDataIdx}.bC,2));
        arrayData{arrayDataIdx}.pairwiseProb.independenceProb = zeros(numel(arrayData),size(arrayData{arrayDataIdx}.bC,1),size(arrayData{arrayDataIdx}.bC,2));
        for comparisonIdx = 1:numel(arrayData)
            normNumStimsWithBothResponse = zeros(size(arrayData{arrayDataIdx}.bC));
            independenceProb = zeros(size(arrayData{arrayDataIdx}.bC));
            for c = 1:size(arrayData{arrayDataIdx}.bE,1)
                for w = 1:1:size(arrayData{arrayDataIdx}.bE,2)
                    % determine window if automatic windowing
                    if(opts.AUTOMATIC_WINDOW) % get high time and low time automatically
                        opts.WINDOW = getAutomaticWindow(arryaData,arrayDataIdx,c,w,opts);
                    end

                    [normNumStimsWithBothResponse(c,w),independenceProb(c,w),difference(c,w)] = computePairwiseProbabilityOfResponse(arrayData,arrayDataIdx,comparisonIdx,c,w,opts);

                end
            end
            
            arrayData{arrayDataIdx}.pairwiseProb.normNumStimsWithBothResponse(comparisonIdx,:,:) = normNumStimsWithBothResponse;
            arrayData{arrayDataIdx}.pairwiseProb.independenceProb(comparisonIdx,:,:) = independenceProb;
            arrayData{arrayDataIdx}.pairwiseProb.difference(comparisonIdx,:,:) = difference;
        end
    end
    
end


function [opts] = configureOpts(optsInput)

    opts.WINDOW = [1,5]/1000;
    opts.AUTOMATIC_WINDOW = 0;
    opts.BASELINE_TIME = 3/1000;
    opts.SUBTRACT_BASELINE = 1;
    opts.COMPUTE_PAIRWISE = 0;
    
    %% check if in opts and optsInput, overwrite if so
    try
        inputFieldnames = fieldnames(optsInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(opts,inputFieldnames{fn}))
               opts.(inputFieldnames{fn}) = optsInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
end

%% this function computes the probability of response for a single unit
function [meanSpikesEvoked_all, normNumStimsResponsive, meanSpikesEvoked_responsive] = computeProbabilityOfResponse(arrayData,arrayDataIdx,c,w,opts)

    % get numSpikesEvoked/numStims for all stimulations
    baseline = sum(arrayData{arrayDataIdx}.spikeTrialTimes{c,w} < min(opts.BASELINE_TIME))/(-1*arrayData{arrayDataIdx}.bE{c,w}(1)/1000+opts.BASELINE_TIME); % mean spikes/ms
    ROI = sum(arrayData{arrayDataIdx}.spikeTrialTimes{c,w} > opts.WINDOW(1) & arrayData{arrayDataIdx}.spikeTrialTimes{c,w} < opts.WINDOW(2))/(opts.WINDOW(2)-opts.WINDOW(1)); % mean spikes/ms
    if(opts.SUBTRACT_BASELINE)
        meanSpikesEvoked_all = (ROI - baseline)*(opts.WINDOW(2)-opts.WINDOW(1))/(arrayData{arrayDataIdx}.numStims(c,w)-baseline*diff(opts.WINDOW));
    else
        meanSpikesEvoked_all = (ROI)*(opts.WINDOW(2)-opts.WINDOW(1))/(arrayData{arrayDataIdx}.numStims(c,w));
    end
    % get numSpikesEvoked/numStims for only those with a
    % response
    % get num stimulations with response / num stimulations
    stimsWithResponse = unique(arrayData{arrayDataIdx}.stimData{c,w}(arrayData{arrayDataIdx}.spikeTrialTimes{c,w} > opts.WINDOW(1) & arrayData{arrayDataIdx}.spikeTrialTimes{c,w} < opts.WINDOW(2)));

    if(opts.SUBTRACT_BASELINE)
        normNumStimsResponsive = (numel(stimsWithResponse)-baseline*diff(opts.WINDOW))/(arrayData{arrayDataIdx}.numStims(c,w)-baseline*diff(opts.WINDOW));
    else
        normNumStimsResponsive = (numel(stimsWithResponse))/(arrayData{arrayDataIdx}.numStims(c,w));
    end

    if(opts.SUBTRACT_BASELINE)
        meanSpikesEvoked_responsive = meanSpikesEvoked_all*(arrayData{arrayDataIdx}.numStims(c,w)-baseline*diff(opts.WINDOW))/(numel(stimsWithResponse)-baseline*diff(opts.WINDOW));
    else
        meanSpikesEvoked_responsive = meanSpikesEvoked_all*arrayData{arrayDataIdx}.numStims(c,w)/(numel(stimsWithResponse));
    end

end

%% this function computes pairwise probabilities
function [normNumStimsWithBothResponse,independenceProb,difference] = computePairwiseProbabilityOfResponse(arrayData,arrayDataIdx,comparisonIdx,c,w,opts)
    if(arrayDataIdx == comparisonIdx)
        normNumStimsWithBothResponse = -1;
        independenceProb = -1;
        difference = -1;
    else
        %% compute P1*P2 and store
        independenceProb = arrayData{arrayDataIdx}.singleProb.normNumStimsResponsive(c,w)*arrayData{comparisonIdx}.singleProb.normNumStimsResponsive(c,w);

        %% find proportion of stimulations with spikes for both units
        stimsWithResponse{1} = unique(arrayData{arrayDataIdx}.stimData{c,w}(arrayData{arrayDataIdx}.spikeTrialTimes{c,w} > opts.WINDOW(1) & arrayData{arrayDataIdx}.spikeTrialTimes{c,w} < opts.WINDOW(2)));
        stimsWithResponse{2} = unique(arrayData{comparisonIdx}.stimData{c,w}(arrayData{comparisonIdx}.spikeTrialTimes{c,w} > opts.WINDOW(1) & arrayData{comparisonIdx}.spikeTrialTimes{c,w} < opts.WINDOW(2)));

        stimsWithBothResponse = intersect(stimsWithResponse{1},stimsWithResponse{2});

        baseline = sum(arrayData{arrayDataIdx}.spikeTrialTimes{c,w} < min(opts.BASELINE_TIME))/(-1*arrayData{arrayDataIdx}.bE{c,w}(1)/1000+opts.BASELINE_TIME);
        baseline = baseline*sum(arrayData{comparisonIdx}.spikeTrialTimes{c,w} < min(opts.BASELINE_TIME))/(-1*arrayData{comparisonIdx}.bE{c,w}(1)/1000+opts.BASELINE_TIME);

        if(opts.SUBTRACT_BASELINE)
            normNumStimsWithBothResponse = (numel(stimsWithBothResponse)-baseline*diff(opts.WINDOW))/(arrayData{arrayDataIdx}.numStims(c,w)-baseline*diff(opts.WINDOW));
        else
            normNumStimsWithBothResponse = (numel(stimsWithBothResponse))/(arrayData{arrayDataIdx}.numStims(c,w));
        end
        difference = normNumStimsWithBothResponse - independenceProb;
    end
end

%% automatic window finder
function [window] = getAutomaticWindow(arrayData,arrayDataIdx,c,w,opts)

    % find bin with edge at zero 
    zeroBin = find(arrayData{arrayDataIdx}.bE{c,w} == 0);

    % find peak in PSTH between 0 and 5ms
    [~,peakIdx] = max(arrayData{arrayDataIdx}.bC{c,w}(zeroBin+floor(opts.WINDOW(1)/binSize):zeroBin+floor(opts.WINDOW(2)/binSize)));
    peakIdx = peakIdx + zeroBin - 1; % shift over so it corresponds with normal indexing

    % get baseline rate
    baselineRate = mean(arrayData{arrayDataIdx}.bC{c,w}(1:zeroBin-floor(opts.BASELINE_TIME/binSize))); % preTime to -3ms

    % set high and low time based on baseline firing rate
    lowTemp = find(arrayData{arrayDataIdx}.bC{c,w}(1:peakIdx) < baselineRate);
    window(1) = arrayData{arrayDataIdx}.bE{c,w}(lowTemp(end));
    highTemp = find(arrayData{arrayDataIdx}.bC{c,w}(peakIdx:end) < baselineRate);
    window(2) = highTemp(1);
    %             lowTime = bE(peakIdx)-0.001;
    %             highTime = bE(peakIdx)+0.001;
    
end