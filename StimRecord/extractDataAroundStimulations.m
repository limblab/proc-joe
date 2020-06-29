function [ arrayData ] = extractDataAroundStimulations( inputData, fileList, stimInfoFileList,opts )
% this function gets the spike time and bin count information around the
% stimulations. opts provides inputs

    %% configure opts and set default values
    opts = configureOpts(opts);
    arrayData = {};
    amplitude_master_list = []; % 
    %% load in mapfile
    MAP_DATA = loadMapFile(inputData.mapFileName(8:end));
    POS_LIST = [MAP_DATA.row,MAP_DATA.col];
    
    %% for each file, extract data and then combine into a single one
    for fileNumber = 1:numel(fileList)
        disp(fileList(fileNumber).name)
        cd(inputData.folderpath)
        cds = commonDataStructure();
        load(stimInfoFileList(fileNumber).name);
        if(~iscell(stimInfo.chanSent))
            stimInfo.chanSent = mat2cell(stimInfo.chanSent',ones(numel(stimInfo.chanSent),1));
            warning('made chan sent a cell array');
        end
        cd(inputData.folderpath)
        cds.file2cds([inputData.folderpath fileList(fileNumber).name],inputData.task,inputData.ranBy,...
            inputData.monkey,inputData.labnum,inputData.array1,inputData.mapFileName); % DO NOT USE RECOVER PRE SYNC, currently this shifts the units and the analog signal differently

        % check to make sure stimInfo.() are all column vectors
        fnames = fieldnames(stimInfo);
        for fIdx = 1:numel(fnames)
            if(size(stimInfo.(fnames{fIdx}),1) < size(stimInfo.(fnames{fIdx}),2))
                stimInfo.(fnames{fIdx}) = stimInfo.(fnames{fIdx})';
            end
        end
        
        opts.NEURON_NUMBER_ALL = find([cds.units.ID]~=0 & [cds.units.ID]~=255);
%         opts.NEURON_NUMBER_ALL = find([cds.units.ID]~=255);
        %% find stim times based on sync line, store in stimInfo
        if(opts.USE_ANALOG_FOR_STIM_TIMES)
            for j=1:numel(cds.analog)
                syncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,opts.ANALOG_SYNC_LINE));
                if ~isempty(syncIdx)
                    aIdx=j;
                end
            end
            stimInfo.stimOn = cds.analog{aIdx}.t(find(diff(cds.analog{aIdx}.(opts.ANALOG_SYNC_LINE)-mean(cds.analog{aIdx}.(opts.ANALOG_SYNC_LINE))>3)>.5));
            
            if(opts.DOWNSAMPLE_STIM_TIMES && opts.STIMULATIONS_PER_TRAIN == 1)
                idx_keep = find(diff(stimInfo.stimOn) > 1);
                idx_keep = [1;idx_keep + 1];
                stimInfo.stimOn = stimInfo.stimOn(idx_keep);
            end
        end

        %% setup useful variables
        % if opts.USE_STIM_CODE, then rebuild
        % stimInfo.waveSent/stimInfo.chanSent
        if(opts.USE_STIM_CODE)
            stim_codes = [cds.trials.stimCode];
            start_time = [cds.trials.startTime];
            stim_info_mask = ones(numel(stimInfo.stimOn),1);
            stimInfo.waveSent = ones(size(stimInfo.stimOn));
            stimInfo.chanSent = ones(size(stimInfo.stimOn));
            
            for st = 1:numel(stimInfo.stimOn)
                % find most recent stim code
                trial_idx = find(stimInfo.stimOn(st) > start_time,1,'last');
                if(isempty(trial_idx))
                    stim_info_mask(st) = 0;
                else
                    stimInfo.waveSent(st) = stim_codes(trial_idx)+1;
                end
            end
            stimInfo.stimOn = stimInfo.stimOn(stim_info_mask==1);
            stimInfo.stimOff = stimInfo.stimOff(stim_info_mask==1);
            stimInfo.chanSent = stimInfo.chanSent(stim_info_mask==1);
            stimInfo.chanSent(:) = opts.STIM_ELECTRODE{:};
            stimInfo.chanSent = mat2cell(stimInfo.chanSent(:),ones(size(stimInfo.chanSent,1),1),ones(size(stimInfo.chanSent,2),1));
            stimInfo.waveSent = stimInfo.waveSent(stim_info_mask==1);
            
            NUM_WAVEFORM_TYPES = 1;
            WAVEFORM_NUMS = 1;
            NUM_CHANS = 1;
            CHAN_LIST = opts.STIM_ELECTRODE; % placeholder because this isn't given that info
        else
            if(~isempty(opts.NUM_WAVEFORM_TYPES))
                NUM_WAVEFORM_TYPES = opts.NUM_WAVEFORM_TYPES;
                WAVEFORM_NUMS = 1:1:NUM_WAVEFORM_TYPES;
            elseif(any(isfield(stimInfo,'waveSent')))
                NUM_WAVEFORM_TYPES = numel(unique(stimInfo.waveSent));
                WAVEFORM_NUMS = unique(stimInfo.waveSent);
            else
                NUM_WAVEFORM_TYPES = 1;
                WAVEFORM_NUMS = 1;
            end

            if(~isempty(opts.CHAN_LIST))
                NUM_CHANS = numel(opts.CHAN_LIST);
                CHAN_LIST = opts.CHAN_LIST;
            elseif(any(isfield(stimInfo,'chanSent')))
                if(iscell(stimInfo.chanSent))
                    NUM_CHANS = numel(uniquecell(stimInfo.chanSent));
                    CHAN_LIST = uniquecell(stimInfo.chanSent);
                else
                    NUM_CHANS = numel(unique(stimInfo.chanSent));
                    CHAN_LIST = unique(stimInfo.chanSent);
                end
            else
                CHAN_LIST = opts.STIM_ELECTRODE;
                NUM_CHANS = 1;
            end
            
            if(iscell(stimInfo.chanSent) && sum(stimInfo.chanSent{1} == -1) == 1) % try to get chan stim from filename, % handles multiple channels
                fname = fileList(fileNumber).name;
                underscoreIdx = strfind(fname,'_');
                if(~isempty(strfind(fname,'chanStim')))
                    strIdx = strfind(fname,'chanStim');
                    underscoreIdx = underscoreIdx(find(underscoreIdx < strIdx,1,'last'));
                    chanStim = str2num(fname(underscoreIdx+1:strIdx-1));
                elseif(~isempty(strfind(fname,'stim')) && ~isempty(strfind(fname,'chan')))
                    strIdx = [strfind(fname,'chan'), strfind(fname,'stim')];
                    chanStim = str2num(fname(strIdx(end-1)+4:strIdx(end)-1)); % some files have chans....chan#stim
                end
                
                if(~isempty(chanStim))
                    for chanSentIdx = 1:numel(stimInfo.chanSent)
                        stimInfo.chanSent{chanSentIdx} = chanStim;
                    end
                end
                
            end
            
            if(opts.FLAG_WAVEFORM) % try to get stim amp from file name, get idx for master list of waveforms
                fname = fileList(fileNumber).name;
                underscoreIdx = strfind(fname,'_');
                strIdx = strfind(fname,'uA');
                underscoreIdx = underscoreIdx(find(underscoreIdx < strIdx,1,'last'));
                ampStim = str2num(fname(underscoreIdx+1:strIdx-1));
                amplitude_master_list_idx = find(amplitude_master_list == ampStim);
                if(isempty(amplitude_master_list_idx))
                    amplitude_master_list(end+1) = ampStim;
                    amplitude_master_list_idx = numel(amplitude_master_list);
                end
                
                stimInfo.waveSent(:) = amplitude_master_list_idx;
            end
        end
    
        for nn_all_idx = 1:numel(opts.NEURON_NUMBER_ALL)
            opts.NEURON_NUMBER = opts.NEURON_NUMBER_ALL(nn_all_idx);
            
            %% find array data idx based on chan and unit #
            arrayDataIdx = [];
            for arr_idx = 1:size(arrayData,2)
                if(cds.units(opts.NEURON_NUMBER).chan == arrayData{arr_idx}.CHAN_REC && ...
                        cds.units(opts.NEURON_NUMBER).ID == arrayData{arr_idx}.ID)
                    arrayDataIdx = arr_idx;
                end
            end
            
            if(isempty(arrayDataIdx))
                arrayDataIdx = size(arrayData,2) + 1;
            end
            
            %% setup arrays
            spikeTrialTimes = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
            spikeTrueTimes = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
            stimuliData = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
            binEdges = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
            binCounts = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
            kinData = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
            forceData = cell(NUM_CHANS,NUM_WAVEFORM_TYPES);
            
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
                    wave_num = WAVEFORM_NUMS(wave);
                    % find artifacts from channel and wave combination
                    stimsPerTrainMask = zeros(min(numel(stimInfo.stimOn),numel(stimInfo.waveSent)),1);
                    stimsPerTrainMask(1:opts.STIMULATIONS_PER_TRAIN:end) = 1;
%                     stimIdx = find(stimsPerTrainMask == 1 & stimInfo.waveSent(1:numel(stimsPerTrainMask)) == wave & stimInfo.chanSent{1:numel(stimsPerTrainMask)}]' == CHAN_LIST(chan));
                    stimIdx = find(stimsPerTrainMask == 1 & stimInfo.waveSent(1:numel(stimsPerTrainMask)) == wave_num & ...
                        checkChanListEquality(stimInfo.chanSent(1:numel(stimsPerTrainMask)),CHAN_LIST{chan}));

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

                        if(arraySize(chan,wave) > 3/4*size(spikeTrialTimes{chan,wave},2))
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
%                     firingRateStimuli = zeros(sum(stimInfo.waveSent == wave & checkChanListEquality(stimInfo.chanSent(1:numel(stimsPerTrainMask)),CHAN_LIST{chan})),2);
%                     for i = 1:sum(stimInfo.waveSent == wave & checkChanListEquality(stimInfo.chanSent(1:numel(stimsPerTrainMask)),CHAN_LIST{chan}))
%                         firingRateStimuli(i,1) = numel(find(stimuliData{chan,wave} == i & spikeTrialTimes{chan,wave} > -opts.PRE_TIME & spikeTrialTimes{chan,wave} < 0))/(opts.PRE_TIME);
%                         firingRateStimuli(i,2) = numel(find(stimuliData{chan,wave} == i & spikeTrialTimes{chan,wave} < 5/1000 & spikeTrialTimes{chan,wave} > 0.5/1000))/(4.5/1000);
%                     end

        %             binCountsVar{chan,wave} = [mean(firingRateStimuli(:,1)),mean(firingRateStimuli(:,2));...
        %                 var(firingRateStimuli(:,1)),var(firingRateStimuli(:,2))];
                end
            end

            %% grab kin data for each stimulation
            if(opts.GET_KIN && nn_all_idx == 1)
                % if cds has a kin entry, use it
                if(isprop(cds,'kin') && ~isempty(cds.kin))
                    kin_dt = mode(diff(cds.kin.t));
                    for chan = 1:NUM_CHANS
                        for wave = 1:NUM_WAVEFORM_TYPES
                            wave_num = WAVEFORM_NUMS(wave);
                            stimIdx = find(stimInfo.waveSent(1:opts.STIMULATIONS_PER_TRAIN:end) == wave_num & ...
                                checkChanListEquality(stimInfo.chanSent(1:opts.STIMULATIONS_PER_TRAIN:numel(stimsPerTrainMask)),CHAN_LIST{chan}));
                            kinData{chan,wave}.x = nan(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                            kinData{chan,wave}.y = nan(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                            kinData{chan,wave}.vx = nan(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                            kinData{chan,wave}.vy = nan(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                            kinData{chan,wave}.ax = nan(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                            kinData{chan,wave}.ay = nan(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));

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
                    
                else
                    kin_dt = 0.001; % assume this and make nan entries to keep indexing correct
                    for chan = 1:NUM_CHANS
                        for wave = 1:NUM_WAVEFORM_TYPES   
                            wave_num = WAVEFORM_NUMS(wave);
                            stimIdx = find(stimInfo.waveSent(1:opts.STIMULATIONS_PER_TRAIN:end) == wave_num & ...
                                checkChanListEquality(stimInfo.chanSent(1:opts.STIMULATIONS_PER_TRAIN:numel(stimsPerTrainMask)),CHAN_LIST{chan}));
                            
                            kinData{chan,wave}.x = nan(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                            kinData{chan,wave}.y = nan(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                            kinData{chan,wave}.vx = nan(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                            kinData{chan,wave}.vy = nan(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                            kinData{chan,wave}.ax = nan(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                            kinData{chan,wave}.ay = nan(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/kin_dt,1));
                        end
                    end
                end
            end

            
            %% grab force data for each stimulation
            if(opts.GET_FORCE && isprop(cds,'force') && ~isempty(cds.force))
                force_dt = mode(diff(cds.force.t));
                for chan = 1:NUM_CHANS
                    for wave = 1:NUM_WAVEFORM_TYPES
                        wave_num = WAVEFORM_NUMS(wave);
                        stimIdx = find(stimInfo.waveSent(1:opts.STIMULATIONS_PER_TRAIN:end) == wave_num & ...
                            checkChanListEquality(stimInfo.chanSent(1:numel(stimsPerTrainMask)),CHAN_LIST{chan}));
                        forceData{chan,wave}.x = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/force_dt,1));
                        forceData{chan,wave}.y = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/force_dt,1));
                        forceData{chan,wave}.vx = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/force_dt,1));
                        forceData{chan,wave}.vy = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/force_dt,1));
                        forceData{chan,wave}.ax = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/force_dt,1));
                        forceData{chan,wave}.ay = zeros(numel(stimIdx),round((opts.POST_TIME + opts.PRE_TIME)/force_dt,1));

                        for st = 1:numel(stimIdx)
                            forceIdx = find(stimInfo.stimOn(stimIdx(st)) <= cds.force.t,1,'first');
                            if(~isempty(forceIdx))
                                forceData{chan,wave}.x(st,:) = cds.kin.x(round(kinIdx-opts.PRE_TIME/force_dt,1):round(kinIdx+opts.POST_TIME/force_dt-1,1));
                                forceData{chan,wave}.y(st,:) = cds.kin.y(round(kinIdx-opts.PRE_TIME/force_dt,1):round(kinIdx+opts.POST_TIME/force_dt-1,1));
                                forceData{chan,wave}.vx(st,:) = cds.kin.vx(round(kinIdx-opts.PRE_TIME/force_dt,1):round(kinIdx+opts.POST_TIME/force_dt-1,1));
                                forceData{chan,wave}.vy(st,:) = cds.kin.vy(round(kinIdx-opts.PRE_TIME/force_dt,1):round(kinIdx+opts.POST_TIME/force_dt-1,1));
                                forceData{chan,wave}.ax(st,:) = cds.kin.ax(round(kinIdx-opts.PRE_TIME/force_dt,1):round(kinIdx+opts.POST_TIME/force_dt-1,1));
                                forceData{chan,wave}.ay(st,:) = cds.kin.ay(round(kinIdx-opts.PRE_TIME/force_dt,1):round(kinIdx+opts.POST_TIME/force_dt-1,1));
                            end
                        end
                    end
                end
            end
            
            if(size(arrayData,2) < arrayDataIdx) % append data
                %% store unit data
                arrayData{arrayDataIdx}.spikeTrialTimes = spikeTrialTimes;
                arrayData{arrayDataIdx}.stimData = stimuliData;
                arrayData{arrayDataIdx}.numStims = numStims;
                arrayData{arrayDataIdx}.binCounts = binCounts;
                arrayData{arrayDataIdx}.binEdges = binEdges;

                if(arrayDataIdx == 1)
                    arrayData{arrayDataIdx}.kin = kinData;
                    arrayData{arrayDataIdx}.force = forceData;
                end
                %% save useful stimulation related information
                arrayData{arrayDataIdx}.CHAN_LIST = CHAN_LIST;
                arrayData{arrayDataIdx}.STIM_PARAMETERS = stimInfo.parameters;
                arrayData{arrayDataIdx}.WAVEFORM_SENT = stimInfo.waveSent;
                arrayData{arrayDataIdx}.CHAN_SENT = stimInfo.chanSent;
                arrayData{arrayDataIdx}.CHAN_REC = cds.units(opts.NEURON_NUMBER).chan;
                arrayData{arrayDataIdx}.WAVEFORM_LIST = WAVEFORM_NUMS;
                
                %% save useful array information
                arrayData{arrayDataIdx}.ID = cds.units(opts.NEURON_NUMBER).ID;
                arrayData{arrayDataIdx}.NN = opts.NEURON_NUMBER;
                
                posIdx = find(MAP_DATA.chan==cds.units(opts.NEURON_NUMBER).chan);
                arrayData{arrayDataIdx}.ROW = 11 - POS_LIST(posIdx,1);
                arrayData{arrayDataIdx}.COL = POS_LIST(posIdx,2);
                
                
            else % merge data with prexisting data
                %% append unit data
                for chan = 1:NUM_CHANS
                    for wave = 1:NUM_WAVEFORM_TYPES
                        arrayData{arrayDataIdx}.spikeTrialTimes{chan,wave} = [arrayData{arrayDataIdx}.spikeTrialTimes{chan,wave},spikeTrialTimes{chan,wave}];
                        arrayData{arrayDataIdx}.stimData{chan,wave} = [arrayData{arrayDataIdx}.stimData{chan,wave},stimuliData{chan,wave}+arrayData{arrayDataIdx}.numStims(chan,wave)];
                        arrayData{arrayDataIdx}.numStims(chan,wave) = arrayData{arrayDataIdx}.numStims(chan,wave)+numStims(chan,wave);
                        arrayData{arrayDataIdx}.binCounts{chan,wave} = arrayData{arrayDataIdx}.binCounts{chan,wave}+binCounts{chan,wave};
                        if(arrayDataIdx == 1)
                            arrayData{arrayDataIdx}.kin{chan,wave}.x = [arrayData{arrayDataIdx}.kin{chan,wave}.x;kinData{chan,wave}.x];
                            arrayData{arrayDataIdx}.kin{chan,wave}.y = [arrayData{arrayDataIdx}.kin{chan,wave}.y;kinData{chan,wave}.y];
                            arrayData{arrayDataIdx}.kin{chan,wave}.vx = [arrayData{arrayDataIdx}.kin{chan,wave}.vx;kinData{chan,wave}.vx];
                            arrayData{arrayDataIdx}.kin{chan,wave}.vy = [arrayData{arrayDataIdx}.kin{chan,wave}.vy;kinData{chan,wave}.vy];
                            arrayData{arrayDataIdx}.kin{chan,wave}.ax = [arrayData{arrayDataIdx}.kin{chan,wave}.ax;kinData{chan,wave}.ax];
                            arrayData{arrayDataIdx}.kin{chan,wave}.ay = [arrayData{arrayDataIdx}.kin{chan,wave}.ay;kinData{chan,wave}.ay];

                            arrayData{arrayDataIdx}.force{chan,wave} = [arrayData{arrayDataIdx}.force{chan,wave};forceData{chan,wave}];
                        end
                    end
                end
                %% append useful stim related info
                arrayData{arrayDataIdx}.WAVEFORM_SENT = [arrayData{arrayDataIdx}.WAVEFORM_SENT;stimInfo.waveSent];
                arrayData{arrayDataIdx}.CHAN_SENT = [arrayData{arrayDataIdx}.CHAN_SENT;stimInfo.chanSent];
            end
            
        end
        
    end

    %% prune waveform/chan combos that have 0 stims
    for arrayDataIdx = 1:numel(arrayData)
        arrayData_mask = ones(size(arrayData{arrayDataIdx}.binCounts));
        for chan = 1:size(arrayData{arrayDataIdx}.binCounts,1)
            for wave = 1:size(arrayData{arrayDataIdx}.binCounts,2)
                if(arrayData{arrayDataIdx}.numStims(chan,wave) == 0)
                    arrayData_mask(chan,wave) = 0;
                end
            end
        end
        if(~sum(sum(arrayData_mask)) == numel(arrayData_mask)) % prune
            arrayData{arrayDataIdx}.spikeTrialTimes = arrayData{arrayDataIdx}.spikeTrialTimes(arrayData_mask==1);
            arrayData{arrayDataIdx}.stimData = arrayData{arrayDataIdx}.stimData(arrayData_mask==1);
            arrayData{arrayDataIdx}.binCounts = arrayData{arrayDataIdx}.binCounts(arrayData_mask==1);
            arrayData{arrayDataIdx}.binEdges = arrayData{arrayDataIdx}.binEdges(arrayData_mask==1);
            arrayData{arrayDataIdx}.kin = arrayData{arrayDataIdx}.kin(arrayData_mask==1);
            arrayData{arrayDataIdx}.force = arrayData{arrayDataIdx}.force(arrayData_mask == 1);
            
            arrayData{arrayDataIdx}.numStims = arrayData{arrayDataIdx}.numStims(arrayData_mask == 1);
            temp = repmat(arrayData{arrayDataIdx}.CHAN_LIST,1,NUM_WAVEFORM_TYPES)
            arrayData{arrayDataIdx}.STIM_PARAM_LIST(:,1) = temp(arrayData_mask == 1);
            if(~isempty(amplitude_master_list))
                temp = repmat(amplitude_master_list,NUM_CHANS,1);
                arrayData{arrayDataIdx}.STIM_PARAM_LIST(:,2) = temp(arrayData_mask == 1);
                arrayData{arrayDataIdx}.CHAN_LIST = arrayData{arrayDataIdx}.CHAN_LIST(find(sum(arrayData_mask,2) > 0));
            end
        end
    end
    
    %% adjust binCounts based on number of stims
    for arrayDataIdx = 1:numel(arrayData)
        arrayData{arrayDataIdx}.binMaxYLim = 0;
        for chan = 1:size(arrayData{arrayDataIdx}.binCounts,1)
            for wave = 1:size(arrayData{arrayDataIdx}.binCounts,2)
                arrayData{arrayDataIdx}.binCounts{chan,wave} = arrayData{arrayDataIdx}.binCounts{chan,wave}/arrayData{arrayDataIdx}.numStims(chan,wave);
                arrayData{arrayDataIdx}.binMaxYLim = max(arrayData{arrayDataIdx}.binMaxYLim,max(arrayData{arrayDataIdx}.binCounts{chan,wave}));
            end
        end
        arrayData{arrayDataIdx}.binMaxYLim = arrayData{arrayDataIdx}.binMaxYLim*1.1;
        
        % put in amp master list that correspond to waveform
        arrayData{arrayDataIdx}.amp_list = amplitude_master_list;

    end
end

function [opts] = configureOpts(optsInput)

    opts.ALIGN_WAVES = 1;
    opts.STIMULI_RESPONSE = 'all'; %'all', 'responsive' or 'nonresponsive'

    opts.STIM_ELECTRODE = -1;
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
    opts.GET_FORCE = 0;
    opts.USE_ANALOG_FOR_STIM_TIMES = 1;
    opts.DOWNSAMPLE_STIM_TIMES = 0;
    opts.ANALOG_SYNC_LINE = 'ainp16';
    opts.USE_STIM_CODE = 0;

    opts.FLAG_WAVEFORM = 0;
    opts.CHAN_LIST = [];
    opts.NUM_WAVEFORM_TYPES = [];
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
