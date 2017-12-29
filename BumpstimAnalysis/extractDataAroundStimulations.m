function [ outputData ] = extractDataAroundStimulations( cds,opts )
% this function gets the spike time and bin count information around the
% stimulations. opts provides inputs

    %% configure opts and set default values
    opts = configureOpts(opts);

    
    %% setup useful variables
    if(any(isfield(cds,'waveforms')))
        numWaveformTypes = numel(unique(cds.waveforms.waveSent));
    else
        numWaveformTypes = 1;
    end

    if(any(isfield(cds.waveforms,'chanSent')))
        numChans = numel(unique(cds.waveforms.chanSent));
        chanList = unique(cds.waveforms.chanSent);
    else
        chanList = opts.STIM_ELECTRODE;
        numChans = 1;
    end
    
    flagResponse = 1;
    
    %% setup arrays
    spikeTrialTimes = cell(numChans,numWaveformTypes);
    spikeTrueTimes = cell(numChans,numWaveformTypes);
    stimuliData = cell(numChans,numWaveformTypes);
    binEdges = cell(numChans,numWaveformTypes);
    binCounts = cell(numChans,numWaveformTypes);
    binCountsVar = cell(numChans,numWaveformTypes);
    
    for c = 1:numChans
        for i = 1:numWaveformTypes
            spikeTrialTimes{c,i} = zeros(1,opts.INITIAL_ARRAY_SIZE);
            spikeTrueTimes{c,i} = zeros(1,opts.INITIAL_ARRAY_SIZE);
            stimuliData{c,i} = zeros(1,opts.INITIAL_ARRAY_SIZE);
        end
    end
    arraySize = zeros(numChans,numWaveformTypes);
    stimNum = zeros(numChans,numWaveformTypes);
    
    %% extract raster data by looking at each stimulation
    
    for st = 1:opts.STIMULATIONS_PER_TRAIN:numel(cds.stimOn)
        spikeMask = cds.units(neuronNumber).spikes.ts > cds.stimOn(st)-opts.PRE_TIME & cds.units(neuronNumber).spikes.ts < cds.stimOn(st)+opts.POST_TIME;
        spikesPlot = (cds.units(neuronNumber).spikes.ts(spikeMask) - cds.stimOn(st));
        mask = spikesPlot > 0 & spikesPlot < min(15/1000,timeAfterStimRawArtifact);
        if(sum(mask) == 0) % no spikes in the range, no response
            flagResponse = 0;
        else
            flagResponse = 1; % spikes in range, response
        end

        numWaves = sum(spikeMask==1);
        
        if(opts.ALIGN_WAVES) % add 4/30 ms to get to negative deflection of spike since positive deflection is what we are aligned on
            spikesPlot = spikesPlot + 4/30000;
        end
        
        % find chan and wave number if they exist (they should)
        chanNumber = 1;
        waveNumber = 1;
        if(any(isfield(cds.waveforms,'chanSent'))) % should always be true, older files might not have this though
            chanNumber = find(unique(cds.waveforms.chanSent)==cds.waveforms.chanSent(st));
        end
        if(any(isfield(cds,'waveforms')))
            waveNumber = cds.waveforms.waveSent(st);
        end
        
        % update stimNum and arrays if plotting all stimuli, if there is a response
        % and we are plotting only those with response, or if there is not
        % a response and that is what we are plotting
        if((plotAllStimuli) || (plotOnlyStimuliWithResponse && flagResponse) || (plotOnlyStimuliWithoutResponse && ~flagResponse))
            stimNum(chanNumber,waveNumber) = stimNum(chanNumber,waveNumber) + 1;

            % update spike times 
            if(~isempty(spikesPlot))
                spikeTrialTimes{chanNumber,waveNumber}(arraySize(chanNumber,waveNumber)+1:arraySize(chanNumber,waveNumber)+numWaves) = spikesPlot';
                spikeTrueTimes{chanNumber,waveNumber}(arraySize(chanNumber,waveNumber)+1:arraySize(chanNumber,waveNumber)+numWaves) = spikesPlot' + cds.stimOn(st);
                stimuliData{chanNumber,waveNumber}(arraySize(chanNumber,waveNumber)+1:arraySize(chanNumber,waveNumber)+numWaves) = stimNum(chanNumber,waveNumber);
                arraySize(chanNumber,waveNumber) = arraySize(chanNumber,waveNumber) + numWaves;
                
                % check if array size is too large, if yes expand arrays
                if(arraySize(chanNumber,waveNumber) > 2/3*size(spikeTrialTimes{chanNumber,waveNumber},2))
                    tempArray = zeros(1,size(spikeTrialTimes{chanNumber,waveNumber},2) + opts.ADDITIONAL_ARRAY_SIZE);
                    tempArray(1:arraySize(chanNumber,waveNumber)) = spikeTrialTimes{chanNumber,waveNumber}(1:arraySize(chanNumber,waveNumber));
                    spikeTrialTimes{chanNumber,waveNumber} = tempArray;
                end
            end
        end
        
    end
    
    % prune arrays
    for c = 1:numChans
        for i = 1:numWaveformTypes
            spikeTrialTimes{c,i} = spikeTrialTimes{c,i}(1:arraySize(c,i));
            spikeTrueTimes{c,i} = spikeTrueTimes{c,i}(1:arraySize(c,i));
            stimuliData{c,i} = stimuliData{c,i}(1:arraySize(c,i));
        end
    end
    
    %% bin data
    for chan = 1:numChans
        for wave = 1:numWaveformTypes
            % bin data
            binEdges{chan,wave} = -opts.PRE_TIME:opts.BIN_SIZE:opts.POST_TIME;
            [binCounts{chan,wave},binEdges{chan,wave}] = histcounts(spikeTrialTimes{chan,wave},binEdges{chan,wave});
            binEdges{chan,wave} = binEdges{chan,wave}*1000;
            binCounts{chan,wave} = binCounts{chan,wave}/stimNum(chan,wave);
            % compute a variance for the data from preTime to -2/1000 and for
            % the data from 1.5/1000 to 5/1000
            firingRateStimuli = zeros(sum(cds.waveforms.waveSent == wave & cds.waveforms.chanSent == chanList(chan)),2);
            for i = 1:sum(cds.waveforms.waveSent == wave & cds.waveforms.chanSent == chanList(chan))
                firingRateStimuli(i,1) = numel(find(stimuliData{chan,wave} == i & spikeTrialTimes{chan,wave} > -opts.PRE_TIME & spikeTrialTimes{chan,wave} < 0))/(opts.PRE_TIME);
                firingRateStimuli(i,2) = numel(find(stimuliData{chan,wave} == i & spikeTrialTimes{chan,wave} < 5/1000 & spikeTrialTimes{chan,wave} > 0.5/1000))/(4.5/1000);
            end

            binCountsVar{chan,wave} = [mean(firingRateStimuli(:,1)),mean(firingRateStimuli(:,2));...
                var(firingRateStimuli(:,1)),var(firingRateStimuli(:,2))];
        end
    end
    
    %% store output data
    outputData.spikeTrialTimes = spikeTrialTimes;
    outputData.spikeTrueTimes = spikeTrueTimes;
    outputData.stimData = stimuliData;
    outputData.numStims = stimNum;
    outputData.bC = binCounts;
    outputData.bE = binEdges;
    outputData.bCVar = binCountsVar;
    
end

function [opts] = configureOpts(optsInput)

    opts.ALIGN_WAVES = 1;
    opts.STIMULI_RESPONSE = 'all'; %'all', 'responsive' or 'nonresponsive'

    opts.STIM_ELECTRODE = 1;
    opts.STIMULATIONS_PER_TRAIN = 1;
    opts.INITIAL_ARRAY_SIZE = 3000;
    opts.ADDITIONAL_ARRAY_SIZE = 500;
    
    opts.PRE_TIME = -5/1000;
    opts.POST_TIME = 30/1000;
    opts.BIN_SIZE = 0.2/1000;
    
    %% check if in optsSave and optsSaveInput, overwrite if so
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