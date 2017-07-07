function [outputFigures, outputData ] = processStimArtifactData(folderpath, inputData )
    %script to load stimulation files and generate perievent plots of 30khz
    %data. Formatted to work with runDataProcessing
    outputFigures = [];
    outputData=[];
    %get list of all files in the folder:
    if ~strcmp(folderpath(end),filesep)
        folderpath=[folderpath,filesep];
    end
    cd(folderpath);
    if(any(isfield(inputData,'fileListExtension')))
        fileList=dir(inputData.fileListExtension);
    else
        fileList=dir('*.nev');
    end
    mapData=loadMapFile(inputData.mapFile(8:end));
    for i=1:size(mapData,1)
        mapData.label{i}=[inputData.array1(6:end),mapData.label{i}];
    end
    %establish trackers for information that is common to all files:
    eList={};
    posList=[];
    chList=[];
    % build filter
    [bFilter,aFilter] = butter(6,500/(30000/2),'high');
    filterParams.b = bFilter;
    filterParams.a = aFilter;
    filterParams.order = 6;
    filterParams.type = 'high';
    thresholdMult = 4;
    % variables to store spike information   
    preOffset = 27;
    postOffset = 20;
    lengthWave = preOffset+postOffset+1;
    numZeros = 200;
    
    % cds and extraseconds for merge purposes
    cds = [];
    nevData = [];
    rawData = [];
    artifactDataTime = 10; % in ms
    artifactData.artifact = zeros(3000,96,artifactDataTime*30000/1000);
    artifactData.t = zeros(3000,1);
    artifactDataIndex = 1;
     
    for i=1:numel(fileList)
        %% load file
        disp(['working on:'])
        disp(fileList(i).name)
        
        cdsTemp=commonDataStructure();
        cdsTemp.file2cds([folderpath,fileList(i).name],inputData.ranBy,inputData.array1,inputData.monkey,inputData.lab,'ignoreJumps',inputData.task,inputData.mapFile,'recoverPreSync');
        
        %% find sync signal in analog data
        useSync=true;
        aIdx=[];
        syncIdx=[];
        if isempty(inputData.useSyncLabel)
        %look for sync under the label sync
            for j=1:numel(cdsTemp.analog)
                syncIdx=find(strcmp(cdsTemp.analog{j}.Properties.VariableNames,'sync'));
                if ~isempty(syncIdx)
                    aIdx=j;
                    syncName='sync';
                end
            end
            %if it wasn't called sync, try for matt's 'StimTrig' label:
            if isempty(syncIdx)
                for j=1:numel(cdsTemp.analog)
                    syncIdx=find(strcmp(cdsTemp.analog{j}.Properties.VariableNames,'StimTrig'));
                    if ~isempty(syncIdx)
                        aIdx=j;
                        syncName='StimTrig';
                    end
                end
            end
            %if we didn't find a sync channel, just look for ainp16
            if isempty(syncIdx)
                for j=1:numel(cdsTemp.analog)
                    syncIdx=find(strcmp(cdsTemp.analog{j}.Properties.VariableNames,'ainp16'));
                    if ~isempty(syncIdx)
                        useSync=false;
                        aIdx=j;
                        syncName='ainp16';
                    end
                end
            end
            if isempty(aIdx)
                error('processStimArtifact:cantFindSync','couldnt find a sync signal')
            end
        else
            useSync=inputData.useSyncLabel;
            if useSync
                %find aIdx:
                for j=1:numel(cdsTemp.analog)
                    syncIdx=find(strcmp(cdsTemp.analog{j}.Properties.VariableNames,'sync'));
                    if ~isempty(syncIdx)
                        aIdx=j;
                        syncName='sync';
                    end
                end
            else
                for j=1:numel(cdsTemp.analog)
                    syncIdx=find(strcmp(cdsTemp.analog{j}.Properties.VariableNames,'ainp16'));
                    if ~isempty(syncIdx)
                        useSync=false;
                        aIdx=j;
                        syncName='ainp16';
                    end
                end
            end
        end
        %% use sync to get stim times:
        artifactDataPre.stimOn=find(diff(cdsTemp.analog{aIdx}.(syncName)-mean(cdsTemp.analog{aIdx}.(syncName))>3)>.5);
        stimOff=find(diff(cdsTemp.analog{aIdx}.(syncName)-mean(cdsTemp.analog{aIdx}.(syncName))<-3)>.5);
        artifactDataPre.stimOn = artifactDataPre.stimOn;
        artifactDataPre.stimOff=nan(size(artifactDataPre.stimOn));
        for j=1:numel(artifactDataPre.stimOn)
            if j<numel(artifactDataPre.stimOn)
                next=artifactDataPre.stimOn(j+1);
            else
                next=numel(cdsTemp.analog{aIdx}.(syncName));
            end
            offIdx=stimOff(find((stimOff>artifactDataPre.stimOn(j)& stimOff<next),1,'first'));
            if ~isempty(offIdx)
                artifactDataPre.stimOff(j)=offIdx;
            end
        end
        
        %% extract data from cdsTemp so that the rest of the code moves quicker hopefully
        cdsTempLFP = cdsTemp.lfp{:,:};
        cdsTempDuration = cdsTemp.meta.duration;

        % merge cdsTemp with cds <- this one is everything but lfp + analog
        if(strcmp(inputData.task,'taskCObump'))
            if(i==1)
                % rewrite everything into cds -- which is really not a
                % commonDataStructure object but whatever
                cds.offsets = [0];
                cds.meta = cdsTemp.meta;
                cds.meta.hasLfp = 0;
                cds.kin = cdsTemp.kin;
                cds.force = cdsTemp.force;
                cds.lfp = {};
                cds.emg = cdsTemp.emg;
                cds.analog = {};
                cds.stimOn = artifactDataPre.stimOn/30000;
                cds.stimOff = artifactDataPre.stimOff/30000;
                cds.triggers = cdsTemp.triggers;
                cds.units = cdsTemp.units;
                cds.trials = cdsTemp.trials;
                cds.aliasList = cdsTemp.aliasList;
                cds.operationLog = cdsTemp.operationLog;
                cds.meta.duration = 0;
            else     
                % port over kin data
                tempKin = cdsTemp.kin;
                tempKin.t = tempKin.t + cds.meta.duration;
                cdsKin = cds.kin;
                cdsKin{end+1:end+size(tempKin,1),:} = tempKin{:,:};
                cds.kin = cdsKin;
                % port over stimOn and stimOff data
                cds.stimOn = [cds.stimOn;artifactDataPre.stimOn/30000 + cds.meta.duration];
                cds.stimOff = [cds.stimOff;artifactDataPre.stimOff/30000 + cds.meta.duration];
                % port over trial data
                tempTrials = cdsTemp.trials;
                tempTrials.number = tempTrials.number + cds.trials.number(end);
                tempTrials.startTime = tempTrials.startTime + cds.meta.duration;
                tempTrials.endTime = tempTrials.endTime + cds.meta.duration;
                tempTrials.tgtOnTime = tempTrials.tgtOnTime + cds.meta.duration;
                tempTrials.goCueTime = tempTrials.goCueTime + cds.meta.duration;
                tempTrials.bumpTime = tempTrials.bumpTime + cds.meta.duration;
                tempTrials.stimTime = tempTrials.stimTime + cds.meta.duration;
                cdsTrials = cds.trials;
                cdsTrials(end+1:end+size(cdsTemp.trials,1),:) = tempTrials(:,:);
                cds.trials = cdsTrials;
                % update meta information
                cds.offsets(end+1,1) = cds.meta.duration;
                cds.meta.dataWindow(2) = cds.meta.dataWindow(2) + cdsTemp.meta.dataWindow(2);
                cds.meta.numTrials = cds.meta.numTrials + cdsTemp.meta.numTrials;
                cds.meta.numReward = cds.meta.numReward + cdsTemp.meta.numReward;
                cds.meta.numAbort = cds.meta.numAbort + cdsTemp.meta.numAbort;
                cds.meta.numFail = cds.meta.numFail + cds.meta.numFail;
                cds.meta.numIncomplete = cds.meta.numIncomplete + cds.meta.numIncomplete;
            end
        else % ignore trial table completely
            if(i==1)
                % rewrite everything into cds -- which is really not a
                % commonDataStructure object but whatever
                cds.offsets = [0];
                cds.meta = cdsTemp.meta;
                cds.meta.hasLfp = 0;
                cds.kin = cdsTemp.kin;
                cds.force = cdsTemp.force;
                cds.lfp = {};
                cds.emg = cdsTemp.emg;
                cds.analog = {};
                cds.stimOn = artifactDataPre.stimOn/30000;
                cds.stimOff = artifactDataPre.stimOff/30000;
                cds.triggers = cdsTemp.triggers;
                cds.units = cdsTemp.units;
                cds.aliasList = cdsTemp.aliasList;
                cds.operationLog = cdsTemp.operationLog;
                cds.meta.duration = 0;
            else     
                % port over kin data
                tempKin = cdsTemp.kin;
                tempKin.t = tempKin.t + cds.meta.duration;
                cdsKin = cds.kin;
                cdsKin{end+1:end+size(tempKin,1),:} = tempKin{:,:};
                cds.kin = cdsKin;
                % port over stimOn and stimOff data
                cds.stimOn = [cds.stimOn;artifactDataPre.stimOn/30000 + cds.meta.duration];
                cds.stimOff = [cds.stimOff;artifactDataPre.stimOff/30000 + cds.meta.duration];                
                % update meta information
                cds.offsets(end+1,1) = cds.meta.duration;
                cds.meta.dataWindow(2) = cds.meta.dataWindow(2) + cdsTemp.meta.dataWindow(2);
            end
        end
        
        clear cdsTemp
        %% use stim times to filter backwards
        % fix stim times this time -- probably wont need to do for all
        % other times because we fixed it during experiment
%         stimOnTemp = [];
%         for stimIdx = 1:numel(artifactDataPre.stimOn)
%             stimOnTemp = [stimOnTemp; artifactDataPre.stimOn(stimIdx) + 300*(0:1:9)'];
%         end
%         artifactDataPre.stimOn = stimOnTemp;
%         artifactDataPre.stimOff = artifactDataPre.stimOn + 2;
        
        % get thresholds for each channel based on non stim data
        thresholdAll = zeros(size(cdsTempLFP,2)-1,1);
        for ch = 2:size(cdsTempLFP,2)
            numPoints = 0;
            for stimuli = 1:numel(artifactDataPre.stimOn)+1
                if(stimuli == 1) % all data before first stim
                    stimData = cdsTempLFP(1:artifactDataPre.stimOn(stimuli),ch);
                elseif(stimuli == numel(artifactDataPre.stimOn) + 1) % all data after last stim artifact
                    stimData = cdsTempLFP(artifactDataPre.stimOn(stimuli-1)+3*30:end,ch);
                else % data before ith stim up to ith-1 
                    stimData = cdsTempLFP(artifactDataPre.stimOn(stimuli-1)+3*30:artifactDataPre.stimOn(stimuli),ch);
                end
                numPoints = numPoints + numel(stimData);
                stimDataTemp = [stimData(:,1);zeros(numZeros,1)];
                stimDataTemp = fliplr(filter(bFilter,aFilter,fliplr(stimDataTemp')))';
                thresholdAll(ch-1) = thresholdAll(ch-1) + sum(stimDataTemp.^2);

            end
            thresholdAll(ch-1) = sqrt(thresholdAll(ch-1)/numPoints);
        end
        
        
        spikeWaves = zeros(10000,lengthWave);
        spikeTimes = zeros(10000,1);
        spikeChan = zeros(10000,1);    
        spikeNum = 1;
        disp('filtering and thresholding')
        for stimIdx = 1:numel(artifactDataPre.stimOn)+1
            stimData = [];
            
            if(stimIdx == 1) % all data before first stim
                stimData = cdsTempLFP(1:artifactDataPre.stimOn(stimIdx),2:end);
            elseif(stimIdx == numel(artifactDataPre.stimOn) + 1) % all data after last stim
                stimData = cdsTempLFP(artifactDataPre.stimOn(stimIdx-1):end,2:end);
                stimData(1:30*1,:) = 0; % blank the first millisecond
            else % data before ith stim up to ith-1 
                stimData = cdsTempLFP(artifactDataPre.stimOn(stimIdx-1):artifactDataPre.stimOn(stimIdx),2:end);
                stimData(1:30*1,:) = 0; % blank the first millisecond
            end
            
            % filter backwards on all channels and threshold
            for ch = 1:size(stimData,2)
                stimDataTemp = [stimData(:,ch);zeros(numZeros,1)];
                
                stimDataTemp = fliplr(filter(bFilter,aFilter,fliplr(stimDataTemp')))';
                stimData(:,ch) = stimDataTemp(1:end-numZeros);
                threshold = thresholdMult*thresholdAll(ch);
%                 threshold = abs(rms(stimData(max(1,numel(stimData(:,ch))-10):end,ch))*thresholdMult);
                thresholdCrossings = find(stimData(:,ch)>threshold);
                
                % remove chains -- find best spot
                idx = 2;
                chain = [1];
                crossingsKeep = [];
                while idx <= numel(thresholdCrossings)
                    if(thresholdCrossings(idx) == thresholdCrossings(idx-1)+1) % store in chain
                        chain = [chain;idx];
                    elseif(~isempty(chain)) % broke a chain, store minidx, update idx, empty chain
                        [~,maxIdx] = max(stimData(thresholdCrossings(chain)));
                        if(isempty(crossingsKeep))
                            crossingsKeep = [thresholdCrossings(maxIdx+chain(1)-1)];
                        else
                            crossingsKeep = [crossingsKeep;thresholdCrossings(maxIdx+chain(1)-1)];
                        end
                        chain = [idx];
                    end
                    idx = idx+1;
                end
                if(numel(thresholdCrossings) > 0)
                    thresholdCrossings = [crossingsKeep;thresholdCrossings(end)];
                end
                
                % remove potential artifacts and too close to beginning
                crossingsMask = ones(numel(thresholdCrossings),1);
                for cross = 1:numel(thresholdCrossings)
                    if(stimData(thresholdCrossings(cross),ch) > 1000 || ...
                            ~(thresholdCrossings(cross)+postOffset <= numel(stimData(:,ch)) && thresholdCrossings(cross)-preOffset > 0))
                        crossingsMask(cross) = 0;
                    end
                end
                thresholdCrossings = thresholdCrossings(crossingsMask(:) == 1);
                
                % go through and weed out ones that are too close to each other
                % prioritize backwards in time
                crossingsMask = ones(numel(thresholdCrossings),1);
                for cross = numel(thresholdCrossings):-1:2
                    if(crossingsMask(cross) == 1) % check time beforehand to see if one is too close
                        crossCheck = cross-1;
                        while crossCheck >= 1 && thresholdCrossings(crossCheck) >= thresholdCrossings(cross) - max(preOffset,postOffset)
                            crossingsMask(crossCheck) = 0;
                            crossCheck = crossCheck-1;
                        end
                    end
                end  
                thresholdCrossings = thresholdCrossings(crossingsMask(:) == 1);

                % store thresholdCrossing data
                for cross = 1:numel(thresholdCrossings)
                    if(stimIdx == 1)
                        spikeTimes(spikeNum) = (thresholdCrossings(cross) - 1)/30000; % this is in seconds
                        
                    else
                        spikeTimes(spikeNum) = (artifactDataPre.stimOn(stimIdx-1) + thresholdCrossings(cross) - 1)/30000; % this is in secondss
                        
                    end
                    spikeChan(spikeNum) = ch;
                    spikeWaves(spikeNum,:) = stimData(thresholdCrossings(cross)-preOffset:thresholdCrossings(cross)+postOffset,ch);
                    spikeNum = spikeNum + 1;

                    if(spikeNum > 0.67*numel(spikeTimes))
                        spikeTimes = [spikeTimes;zeros(1000,1)];
                        spikeWaves = [spikeWaves; zeros(1000,lengthWave)];
                        spikeChan = [spikeChan; zeros(1000,1)];
                    end
                end
                 
                
            end % end channel for
        end % end stimOn for
        
        % store spike data
        disp('storing data')
        spikeTimes = spikeTimes(1:spikeNum-1,1);
        spikeWaves = spikeWaves(1:spikeNum-1,:);
        spikeChan = spikeChan(1:spikeNum-1,:);
        [spikeTimes,sortOrder] = sort(spikeTimes);
        spikeWaves = spikeWaves(sortOrder,:);
        spikeChan = spikeChan(sortOrder,:);
        if(i==1)
            nevData.ts = spikeTimes;
            nevData.waveforms = spikeWaves(:,:);
            nevData.elec = spikeChan(:,:);
        else
            nevData.ts(end+1:end+numel(spikeTimes),:) = spikeTimes + cds.meta.duration;
            nevData.waveforms(end+1:end+numel(spikeTimes),:) = spikeWaves(:,:);
            nevData.elec(end+1:end+numel(spikeTimes),:) = spikeChan(:,:);
        end
        % store raw data where a threshold crossing occurs (index of
        % spikeTimes in this file
        if(i==1)
            rawData.ts = spikeTimes;
            rawIdx = ceil(spikeTimes*30000);
            rawData.waveforms = zeros(numel(spikeTimes),preOffset*2+postOffset*2+1);
            for r = 1:numel(spikeTimes)
                if(rawIdx(r)-preOffset*2 > 0 && rawIdx(r)+postOffset*2 <= size(cdsTempLFP,1))
                    rawData.waveforms(r,:) = cdsTempLFP(rawIdx(r)-preOffset*2:rawIdx(r)+postOffset*2,spikeChan(r)+1);
                end
            end
            rawData.elec = spikeChan;
        else
            rawData.ts(end+1:end+numel(spikeTimes),:) = spikeTimes + cds.meta.duration;
            rawIdx = ceil(spikeTimes*30000);
            rOffset = size(rawData.waveforms,1);
            rawData.waveforms(end+1:end+numel(spikeTimes),:) = zeros(numel(spikeTimes),preOffset*2+postOffset*2+1);
            for r = 1:numel(spikeTimes)
                if(rawIdx(r)-preOffset*2 > 0 && rawIdx(r)+postOffset*2 <= size(cdsTempLFP,1))  
                    rawData.waveforms(r+rOffset,:) = cdsTempLFP(rawIdx(r)-preOffset*2:rawIdx(r)+postOffset*2,spikeChan(r)+1);
                end
            end
            rawData.elec(end+1:end+numel(spikeTimes),:) = spikeChan;
        end

        clear spikeTimes
        clear spikeChan
        clear spikeWaves
        
        % store artifact data
        for art = 1:numel(artifactDataPre.stimOn)
            if(artifactDataPre.stimOn(art) + artifactDataTime*30000/1000 <= size(cdsTempLFP,1))
                artifactData.artifact(artifactDataIndex,:,:) = cdsTempLFP(artifactDataPre.stimOn(art):artifactDataPre.stimOn(art)+floor(artifactDataTime*30000/1000)-1,2:end)';
                if(i==1)
                    artifactData.t(artifactDataIndex,1) = artifactDataPre.stimOn(art)/30000;
                else
                    artifactData.t(artifactDataIndex,1) = artifactDataPre.stimOn(art)/30000 + cds.meta.duration;
                end
                artifactDataIndex = artifactDataIndex + 1;
            end
            if(artifactDataIndex >= size(artifactData.artifact,1))
                artifactData.artifact(end+1:end+1000,:,:) = 0;
                artifactData.t(end+1:end+1000,1) = 0;
            end
        end
        
        cds.meta.duration = cds.meta.duration + cdsTempDuration;
       
    end
    cds.artifactData = artifactData;
    cds.rawData = rawData;
    if(any(isfield(inputData,'stimsPerBump')))
        cds.stimsPerBump = inputData.stimsPerBump;
    else
        cds.stimsPerBump = -1; % usually means no bump, can also mean that the user forgot to input stimsPerBump
    end
    cds.filter = filterParams;
    cds.processed = 0;
    save(strcat(fileList(1).name(1:end-4),'_cds.mat'),'cds','-v7.3');
    save(strcat(fileList(1).name(1:end-4),'_nevData.mat'),'nevData','-v7.3');
end

