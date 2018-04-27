function [ arrayData ] = extractDataAroundStimulations( inputData, fileList, stimInfoFileList,opts )
% this function gets the spike time and bin count information around the
% stimulations. opts provides inputs

    %% configure opts and set default values
    opts = configureOpts(opts);
    arrayData = {};
    
    %% load in mapfile
    MAP_DATA = loadMapFile(inputData.mapFileName(8:end));
    POS_LIST = [MAP_DATA.row,MAP_DATA.col];
    
    %% for each file, extract data and then combine into a single one
    for fileNumber = 1:numel(fileList)
        disp(fileList(fileNumber).name)
        cd(inputData.folderpath)
        cds = commonDataStructure();
        load(stimInfoFileList(fileNumber).name);
        cd(inputData.folderpath)
        cds.file2cds([inputData.folderpath fileList(fileNumber).name],inputData.task,inputData.ranBy,...
            inputData.monkey,inputData.labnum,inputData.array1,inputData.mapFileName); % DO NOT USE RECOVER PRE SYNC, currently this shifts the units and the analog signal differently

%         stimInfo.chanSent = outputData.waveforms.chanSent;
%         stimInfo.waveSent = outputData.waveforms.waveSent;
%         stimInfo.parameters = outputData.waveforms.parameters;
%         
        if(isempty(opts.NEURON_NUMBER_ALL))
            opts.NEURON_NUMBER_ALL = find([cds.units.ID]~=0 & [cds.units.ID]~=255);
        end
        
        %% find stim times based on sync line, store in stimInfo
        if(opts.USE_ANALOG_FOR_STIM_TIMES)
            for j=1:numel(cds.analog)
                syncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,opts.ANALOG_SYNC_LINE));
                if ~isempty(syncIdx)
                    aIdx=j;
                end
            end
            stimInfo.stimOn = cds.analog{aIdx}.t(find(diff(cds.analog{aIdx}.(opts.ANALOG_SYNC_LINE)-mean(cds.analog{aIdx}.(opts.ANALOG_SYNC_LINE))>3)>.5));
        end

        %% setup useful variables
        if(any(isfield(stimInfo,'waveSent')))
            NUM_WAVEFORM_TYPES = numel(unique(stimInfo.waveSent));
        else
            NUM_WAVEFORM_TYPES = 1;
        end

        if(any(isfield(stimInfo,'chanSent')))
            NUM_CHANS = numel(unique(stimInfo.chanSent));
            CHAN_LIST = unique(stimInfo.chanSent);
        else
            CHAN_LIST = opts.STIM_ELECTRODE;
            NUM_CHANS = 1;
        end

        flagResponse = 1;

        for arrayDataIdx = 1:numel(opts.NEURON_NUMBER_ALL)
            opts.NEURON_NUMBER = opts.NEURON_NUMBER_ALL(arrayDataIdx);
            
            %% setup arrays
            spikeTrialTimes = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
            spikeTrueTimes = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
            stimuliData = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
            binEdges = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
            binCounts = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
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
                    stimsPerTrainMask = zeros(min(numel(stimInfo.stimOn),numel(stimInfo.waveSent)),1);
                    stimsPerTrainMask(1:opts.STIMULATIONS_PER_TRAIN:end) = 1;
                    stimIdx = find(stimsPerTrainMask == 1& stimInfo.waveSent(1:numel(stimsPerTrainMask)) == wave & stimInfo.chanSent(1:numel(stimsPerTrainMask)) == CHAN_LIST(chan));

                    numStims(chan,wave) = numel(stimIdx);
                    % process stims in batches to improve speed but control memory
                    for st = 1:opts.STIMULATION_BATCH_SIZE:numel(stimIdx)

                        spikeDataTemp = repmat(cds.units(opts.NEURON_NUMBER).spikes.ts,1,min(numel(stimIdx)-st+1,opts.STIMULATION_BATCH_SIZE)) - ...
                            stimInfo.stimOn(stimIdx(st:min(numel(stimIdx),st+opts.STIMULATION_BATCH_SIZE-1)))';
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
            for chan = 1:NUM_CHANS
                for wave = 1:NUM_WAVEFORM_TYPES
                    % bin data
                    binEdges{chan,wave} = -opts.PRE_TIME:opts.BIN_SIZE:opts.POST_TIME;
                    [binCounts{chan,wave},binEdges{chan,wave}] = histcounts(spikeTrialTimes{chan,wave},binEdges{chan,wave});
                    binEdges{chan,wave} = binEdges{chan,wave}*1000;
                    % compute a variance for the data from preTime to -2/1000 and for
                    % the data from 1.5/1000 to 5/1000
                    firingRateStimuli = zeros(sum(stimInfo.waveSent == wave & stimInfo.chanSent == CHAN_LIST(chan)),2);
                    for i = 1:sum(stimInfo.waveSent == wave & stimInfo.chanSent == CHAN_LIST(chan))
                        firingRateStimuli(i,1) = numel(find(stimuliData{chan,wave} == i & spikeTrialTimes{chan,wave} > -opts.PRE_TIME & spikeTrialTimes{chan,wave} < 0))/(opts.PRE_TIME);
                        firingRateStimuli(i,2) = numel(find(stimuliData{chan,wave} == i & spikeTrialTimes{chan,wave} < 5/1000 & spikeTrialTimes{chan,wave} > 0.5/1000))/(4.5/1000);
                    end

        %             binCountsVar{chan,wave} = [mean(firingRateStimuli(:,1)),mean(firingRateStimuli(:,2));...
        %                 var(firingRateStimuli(:,1)),var(firingRateStimuli(:,2))];
                end
            end

            %% grab kin data for each stimulation
            if(opts.GET_KIN && isfield(cds,'kin') && ~isempty(cds.kin))
                kin_dt = mode(diff(cds.kin.t));
                for chan = 1:NUM_CHANS
                    for wave = 1:NUM_WAVEFORM_TYPES
                        stimIdx = find(stimInfo.waveSent(1:opts.STIMULATIONS_PER_TRAIN:end) == wave & stimInfo.chanSent(1:opts.STIMULATIONS_PER_TRAIN:end) == CHAN_LIST(chan));
                        kinData{chan,wave}.x = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                        kinData{chan,wave}.y = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                        kinData{chan,wave}.vx = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                        kinData{chan,wave}.vy = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                        kinData{chan,wave}.ax = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                        kinData{chan,wave}.ay = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));

                        for st = 1:numel(stimIdx)
                            kinIdx = find(stimInfo.stimOn(stimIdx(st)) <= cds.kin.t,1,'first');
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


            if(fileNumber == 1)
                %% store unit data
                arrayData{arrayDataIdx}.spikeTrialTimes = spikeTrialTimes;
                arrayData{arrayDataIdx}.stimData = stimuliData;
                arrayData{arrayDataIdx}.numStims = numStims;
                arrayData{arrayDataIdx}.bC = binCounts;
                arrayData{arrayDataIdx}.bE = binEdges;

                arrayData{arrayDataIdx}.kin = kinData;

                %% save useful stimulation related information
                arrayData{arrayDataIdx}.CHAN_LIST = CHAN_LIST;
                arrayData{arrayDataIdx}.STIM_PARAMETERS = stimInfo.parameters;
                arrayData{arrayDataIdx}.WAVEFORM_SENT = stimInfo.waveSent;
                arrayData{arrayDataIdx}.CHAN_SENT = stimInfo.chanSent;
                arrayData{arrayDataIdx}.CHAN_REC = cds.units(opts.NEURON_NUMBER).chan;
                
                %% save useful array information
                arrayData{arrayDataIdx}.ID = cds.units(opts.NEURON_NUMBER).ID;
                arrayData{arrayDataIdx}.NN = opts.NEURON_NUMBER;
                
                posIdx = find(MAP_DATA.chan==cds.units(opts.NEURON_NUMBER).chan);
                arrayData{arrayDataIdx}.ROW = 11 - POS_LIST(posIdx,1);
                arrayData{arrayDataIdx}.COL = POS_LIST(posIdx,2);
            else
                %% append unit data
                for chan = 1:NUM_CHANS
                    for wave = 1:NUM_WAVEFORM_TYPES
                        arrayData{arrayDataIdx}.spikeTrialTimes{chan,wave} = [arrayData{arrayDataIdx}.spikeTrialTimes{chan,wave},spikeTrialTimes{chan,wave}];
                        arrayData{arrayDataIdx}.stimData{chan,wave} = [arrayData{arrayDataIdx}.stimData{chan,wave},stimuliData{chan,wave}+arrayData{arrayDataIdx}.numStims(chan,wave)];
                        arrayData{arrayDataIdx}.numStims(chan,wave) = arrayData{arrayDataIdx}.numStims(chan,wave)+numStims(chan,wave);
                        arrayData{arrayDataIdx}.bC{chan,wave} = arrayData{arrayDataIdx}.bC{chan,wave}+binCounts{chan,wave};

                        arrayData{arrayDataIdx}.kin{chan,wave} = [arrayData{arrayDataIdx}.kin{chan,wave};kinData];
                        
                    end
                end
                %% append useful stim related info
                arrayData{arrayDataIdx}.WAVEFORM_SENT = [arrayData{arrayDataIdx}.WAVEFORM_SENT;stimInfo.waveSent];
                arrayData{arrayDataIdx}.CHAN_SENT = [arrayData{arrayDataIdx}.CHAN_SENT;stimInfo.chanSent];
            end
            
        end
        
    end

    for arrayDataIdx = 1:numel(arrayData)
        arrayData{arrayDataIdx}.binMaxYLim = 0;
        for chan = 1:NUM_CHANS
            for wave = 1:NUM_WAVEFORM_TYPES
                arrayData{arrayDataIdx}.bC{chan,wave} = arrayData{arrayDataIdx}.bC{chan,wave}/arrayData{arrayDataIdx}.numStims(chan,wave);
                arrayData{arrayDataIdx}.binMaxYLim = max(arrayData{arrayDataIdx}.binMaxYLim,max(arrayData{arrayDataIdx}.bC{chan,wave}));

            end
        end
        arrayData{arrayDataIdx}.binMaxYLim = arrayData{arrayDataIdx}.binMaxYLim*1.1;
    end
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
    
    opts.NEURON_NUMBER_ALL = [];
    
    opts.GET_KIN = 0;
    opts.USE_ANALOG_FOR_STIM_TIMES = 0;
    opts.ANALOG_SYNC_LINE = 'ainp16';
    
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
    
    opts.NEURON_NUMBER = [];
end
