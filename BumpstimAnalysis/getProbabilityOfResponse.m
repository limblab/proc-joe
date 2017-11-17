function [probOut] = getProbabilityOfResponse(arrayData,opts)

    opts = configureOpts(opts);
    binSize = mean(diff(arrayData{1}.bE{1,1})); % in ms
    %% calculate probability of detection
    for arrayDataIdx = 1:numel(arrayData)
        meanSpikesEvoked_all = zeros(size(arrayData{arrayDataIdx}.bC));
        normNumStimsResponsive = zeros(size(arrayData{arrayDataIdx}.bC));
        meanSpikesEvoked_responsive = zeros(size(arrayData{arrayDataIdx}.bC));
        
        
        for c = 1:size(arrayData{arrayDataIdx},1)
            for w = 1:1:size(arrayData{arrayDataIdx},2)
                if(opts.AUTOMATIC_WINDOW) % get high time and low time automatically
                    % find bin with edge at zero 
                    zeroBin = find(arrayData{arrayDataIdx}.bE{c,w} == 0);

                    % find peak in PSTH between 0 and 5ms
                    [~,peakIdx] = max(arrayData{arrayDataIdx}.bC{c,w}(zeroBin+floor(opts.WINDOW(1)/binSize):zeroBin+floor(opts.WINDOW(2)/binSize)));
                    peakIdx = peakIdx + zeroBin - 1; % shift over so it corresponds with normal indexing
                    
                    % get baseline rate
                    baselineRate = arrayData{arrayDataIdx}.bC{c,w}(1:zeroBin-floor(opts.BASELINE_TIME/binSize)); % preTime to -3ms
                    
                    % set high and low time based on baseline firing rate
                    lowTemp = find(arrayData{arrayDataIdx}.bC{c,w}(1:peakIdx) < baselineRate);
                    opts.WINDOW(1) = arrayData{arrayDataIdx}.bE{c,w}(lowTemp(end));
                    highTemp = find(arrayData{arrayDataIdx}.bC{c,w}(peakIdx:end) < baselineRate);
                    opts.WINDOW(2) = highTemp(1);
        %             lowTime = bE(peakIdx)-0.001;
        %             highTime = bE(peakIdx)+0.001;
                end
                
                % get numSpikesEvoked/numStims for all stimulations
                baseline = sum(arrayData{arrayDataIdx}.spikeTrialTimes{c,w} < min(opts.BASELINE_TIME))/(arrayData{arrayDataIdx}.bE{c,w}(1)-opts.BASELINE_TIME); % mean spikes/ms
                ROI = sum(arrayData{arrayDataIdx}.spikeTrialTimes{c,w} > opts.WINDOW(1) & arrayData{arrayDataIdx}.spikeTrialTimes{c,w} < opts.WINDOW(2))/(opts.WINDOW(2)-opts.WINDOW(1)); % mean spikes/ms
                meanSpikesEvoked_all(c,w) = (ROI - baseline)*(opts.WINDOW(2)-opts.WINDOW(1));

                % get numSpikesEvoked/numStims for only those with a
                % response
                % get num stimulations with response / num stimulations
                stimsWithResponse = unique(arrayData{arrayDataIdx}.stimData{c,w}(arrayData{arrayDataIdx}.spikeTrialTimes{c,w} > opts.WINDOW(1) & arrayData{arrayDataIdx}.spikeTrialTimes{c,w} < opts.WINDOW(2)));
                
                normNumStimsResponsive(c,w) = (numel(stimsWithResponse)-baseline*diff(opts.WINDOW))/(arrayData{arrayDataIdx}.numStims(c,w)-baseline*diff(opts.WINDOW));
                
                meanSpikesEvoked_responsive(c,w) = meanSpikesEvoked_all(c,w)*arrayData{arrayDataIdx}.numStims(c,w)/numel(stimsWithResponse);
                
            end
        end
    end
    probOut.meanSpikesEvoked_all = meanSpikesEvoked_all;
    probOut.normNumStimsResponsive = normNumStimsResponsive;
    probOut.meanSpikesEvokedResponsive = meanSpikesEvoked_responsive;
end


function [opts] = configureOpts(optsInput)

    opts.WINDOW = [1,5]/1000;
    opts.AUTOMATIC_WINDOW = 0;
    opts.BASELINE_TIME = 3/1000;
    
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