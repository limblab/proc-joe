function [ outputData ] = extractDataAroundStimulations( cds,opts )
% this function gets the spike time and bin count information around the
% stimulations. opts provides inputs

    %% configure opts and set default values
    opts = configureOpts(opts);
    
    %% setup useful variables
    if(any(isfield(cds,'waveforms')))
        NUM_WAVEFORM_TYPES = numel(unique(cds.waveforms.waveSent));
    else
        NUM_WAVEFORM_TYPES = 1;
    end

    if(any(isfield(cds.waveforms,'chanSent')))
        NUM_CHANS = numel(unique(cds.waveforms.chanSent));
        CHAN_LIST = unique(cds.waveforms.chanSent);
    else
        CHAN_LIST = opts.STIM_ELECTRODE;
        NUM_CHANS = 1;
    end
    
    flagResponse = 1;
    
    %% setup arrays
    spikeTrialTimes = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
    spikeTrueTimes = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
    stimuliData = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
    binEdges = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
    binCounts = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
    binCountsVar = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
    kinData = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
    
    for c = 1:NUM_CHANS
        for i = 1:NUM_WAVEFORM_TYPES
            spikeTrialTimes{c,i} = zeros(1,opts.INITIAL_ARRAY_SIZE);
            spikeTrueTimes{c,i} = zeros(1,opts.INITIAL_ARRAY_SIZE);
            stimuliData{c,i} = zeros(1,opts.INITIAL_ARRAY_SIZE);
        end
    end
    arraySize = zeros(NUM_CHANS,NUM_WAVEFORM_TYPES);
    numStims = zeros(NUM_CHANS,NUM_WAVEFORM_TYPES);
        
    %% extract raster data by looking at each stimulation
    for chan = 1:NUM_CHANS
        for wave = 1:NUM_WAVEFORM_TYPES
            % find artifacts from channel and wave combination
            stimsPerTrainMask = zeros(min(numel(cds.stimOn),numel(cds.waveforms.waveSent)),1);
            stimsPerTrainMask(1:opts.STIMULATIONS_PER_TRAIN:end) = 1;
            stimIdx = find(stimsPerTrainMask == 1& cds.waveforms.waveSent(1:numel(stimsPerTrainMask)) == wave & cds.waveforms.chanSent(1:numel(stimsPerTrainMask)) == CHAN_LIST(chan));
            
            numStims(chan,wave) = numel(stimIdx);
            % process stims in batches to improve speed but control memory
            for st = 1:opts.STIMULATION_BATCH_SIZE:numel(stimIdx)
                
                spikeDataTemp = repmat(cds.units(opts.NEURON_NUMBER).spikes.ts,1,min(numel(stimIdx)-st+1,opts.STIMULATION_BATCH_SIZE)) - ...
                    cds.stimOn(stimIdx(st:min(numel(stimIdx),st+opts.STIMULATION_BATCH_SIZE-1)))';
                spikeDataTrueTemp = repmat(cds.units(opts.NEURON_NUMBER).spikes.ts,1,min(numel(stimIdx)-st+1,opts.STIMULATION_BATCH_SIZE));
                
                spikeMask = spikeDataTemp > -opts.PRE_TIME & spikeDataTemp < opts.POST_TIME;
                
                spikeTimes = spikeDataTemp(spikeMask);% spike times
                spikeDataTrue = spikeDataTrueTemp(spikeMask); % true spike times
                [~,stimNums] = find(spikeMask); % stim num (column data)
                                
                if(opts.ALIGN_WAVES) % add 4/30 ms to get to negative deflection of spike since positive deflection is what we are aligned on
                    spikeTimes = spikeTimes + 4/30000;
                end
                
                
                numWaves = numel(spikeTimes);
                if(numWaves > 0)
                    spikeTrialTimes{chan,wave}(arraySize(chan,wave)+1:arraySize(chan,wave)+numWaves) = spikeTimes';
                    spikeTrueTimes{chan,wave}(arraySize(chan,wave)+1:arraySize(chan,wave)+numWaves) = spikeDataTrue'; % not correct, fix later
                    stimuliData{chan,wave}(arraySize(chan,wave)+1:arraySize(chan,wave)+numWaves) = (stimNums+st-1)';
                    arraySize(chan,wave) = arraySize(chan,wave) + numWaves;
                end
                
                if(arraySize(chan,wave) > 2/3*size(spikeTrialTimes{chan,wave},2))
                    spikeTrialTimes{chan,wave} = [spikeTrialTimes{chan,wave},zeros(1,opts.ADDITIONAL_ARRAY_SIZE)];
                    spikeTrueTimes{chan,wave} = [spikeTrueTimes{chan,wave},zeros(1,opts.ADDITIONAL_ARRAY_SIZE)];
                    stimuliData{chan,wave} = [stimuliData{chan,wave},zeros(1,opts.ADDITIONAL_ARRAY_SIZE)];
                end
                
            end
        end
    end
    
    % prune arrays
    for c = 1:NUM_CHANS
        for i = 1:NUM_WAVEFORM_TYPES
            spikeTrialTimes{c,i} = spikeTrialTimes{c,i}(1:arraySize(c,i));
            spikeTrueTimes{c,i} = spikeTrueTimes{c,i}(1:arraySize(c,i));
            stimuliData{c,i} = stimuliData{c,i}(1:arraySize(c,i));
        end
    end
    
    %% bin data
    maxYLim = 0;
    for chan = 1:NUM_CHANS
        for wave = 1:NUM_WAVEFORM_TYPES
            % bin data
            binEdges{chan,wave} = -opts.PRE_TIME:opts.BIN_SIZE:opts.POST_TIME;
            [binCounts{chan,wave},binEdges{chan,wave}] = histcounts(spikeTrialTimes{chan,wave},binEdges{chan,wave});
            binEdges{chan,wave} = binEdges{chan,wave}*1000;
            binCounts{chan,wave} = binCounts{chan,wave}/numStims(chan,wave);
            maxYLim = max(max(binCounts{chan,wave})*1.1,maxYLim); 
            % compute a variance for the data from preTime to -2/1000 and for
            % the data from 1.5/1000 to 5/1000
            firingRateStimuli = zeros(sum(cds.waveforms.waveSent == wave & cds.waveforms.chanSent == CHAN_LIST(chan)),2);
            for i = 1:sum(cds.waveforms.waveSent == wave & cds.waveforms.chanSent == CHAN_LIST(chan))
                firingRateStimuli(i,1) = numel(find(stimuliData{chan,wave} == i & spikeTrialTimes{chan,wave} > -opts.PRE_TIME & spikeTrialTimes{chan,wave} < 0))/(opts.PRE_TIME);
                firingRateStimuli(i,2) = numel(find(stimuliData{chan,wave} == i & spikeTrialTimes{chan,wave} < 5/1000 & spikeTrialTimes{chan,wave} > 0.5/1000))/(4.5/1000);
            end

            binCountsVar{chan,wave} = [mean(firingRateStimuli(:,1)),mean(firingRateStimuli(:,2));...
                var(firingRateStimuli(:,1)),var(firingRateStimuli(:,2))];
        end
    end
    
    %% grab kin data for each stimulation
    if(opts.GET_KIN && isfield(cds,'kin') && ~isempty(cds.kin))
        kin_dt = mode(diff(cds.kin.t));
        for chan = 1:NUM_CHANS
            for wave = 1:NUM_WAVEFORM_TYPES
                stimIdx = find(cds.waveforms.waveSent(1:opts.STIMULATIONS_PER_TRAIN:end) == wave & cds.waveforms.chanSent(1:opts.STIMULATIONS_PER_TRAIN:end) == CHAN_LIST(chan));
                kinData{chan,wave}.x = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                kinData{chan,wave}.y = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                kinData{chan,wave}.vx = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                kinData{chan,wave}.vy = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                kinData{chan,wave}.ax = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                kinData{chan,wave}.ay = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));

                for st = 1:numel(stimIdx)
                    kinIdx = find(cds.stimOn(stimIdx(st)) <= cds.kin.t,1,'first');
                    if(~isempty(kinIdx))
                        kinData{chan,wave}.x(st,:) = cds.kin.x(round(kinIdx-opts.PRE_TIME/kin_dt,1):round(kinIdx+opts.POST_TIME/kin_dt-1,1));
                        kinData{chan,wave}.y(st,:) = cds.kin.y(round(kinIdx-opts.PRE_TIME/kin_dt,1):round(kinIdx+opts.POST_TIME/kin_dt-1,1));
                        kinData{chan,wave}.vx(st,:) = cds.kin.vx(round(kinIdx-opts.PRE_TIME/kin_dt,1):round(kinIdx+opts.POST_TIME/kin_dt-1,1));
                        kinData{chan,wave}.vy(st,:) = cds.kin.vy(round(kinIdx-opts.PRE_TIME/kin_dt,1):round(kinIdx+opts.POST_TIME/kin_dt-1,1));
                        kinData{chan,wave}.ax(st,:) = cds.kin.ax(round(kinIdx-opts.PRE_TIME/kin_dt,1):round(kinIdx+opts.POST_TIME/kin_dt-1,1));
                        kinData{chan,wave}.ay(st,:) = cds.kin.ay(round(kinIdx-opts.PRE_TIME/kin_dt,1):round(kinIdx+opts.POST_TIME/kin_dt-1,1));
                    end
                end
            end
        end
    end
    
    %% store output data
    outputData.spikeTrialTimes = spikeTrialTimes;
    outputData.spikeTrueTimes = spikeTrueTimes;
    outputData.stimData = stimuliData;
    outputData.numStims = numStims;
    outputData.bC = binCounts;
    outputData.bE = binEdges;
    outputData.bCVar = binCountsVar;
    outputData.binMaxYLim = maxYLim;
    outputData.kin = kinData;
    
    %% save useful stimulation related information
    outputData.CHAN_LIST = CHAN_LIST;
    outputData.STIM_PARAMETERS = cds.waveforms.parameters;
    outputData.WAVEFORM_SENT = cds.waveforms.waveSent;
    outputData.CHAN_SENT = cds.waveforms.chanSent;
    outputData.CHAN_REC = cds.units(opts.NEURON_NUMBER).chan;
end

function [opts] = configureOpts(optsInput)

    opts.ALIGN_WAVES = 1;
    opts.STIMULI_RESPONSE = 'all'; %'all', 'responsive' or 'nonresponsive'

    opts.STIM_ELECTRODE = 1;
    opts.STIMULATIONS_PER_TRAIN = 1;
    opts.INITIAL_ARRAY_SIZE = 6000;
    opts.ADDITIONAL_ARRAY_SIZE = ceil(opts.INITIAL_ARRAY_SIZE*1/3);
    opts.STIMULATION_BATCH_SIZE = 2000;
    
    opts.PRE_TIME = -5/1000;
    opts.POST_TIME = 30/1000;
    opts.BIN_SIZE = 0.2/1000;
    
    opts.TIME_AFTER_STIMULATION_WAVEFORMS = 10/1000;
    
    opts.NEURON_NUMBER = 1;
    
    opts.GET_KIN = 0;
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
