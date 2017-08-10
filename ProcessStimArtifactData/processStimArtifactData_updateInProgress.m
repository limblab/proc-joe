function [outputFigures, outputData ] = processStimArtifactData(folderpath, inputData )
    %script to load stimulation files and generate perievent plots of 30khz
    %data. Formatted to work with runDataProcessing
    outputFigures = [];
    outputData=[];
    %get list of all files in the folder:
    templateSize = inputData.templateSize*30000;
    noSync = 0;
    noSyncIntended = inputData.noSyncIntended;
    
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

    % variables to store spike information   
    preOffset = 27;
    postOffset = 20;
    maxAmplitude = 700; %uV
    lengthWave = preOffset+postOffset+1;
    numZeros = 200;
    
    % cds and extraseconds for merge purposes
    cds = [];
    nevData = [];
    rawData = [];
    artifactDataTime = inputData.artifactDataTime; % in ms
    artifactData.artifact = zeros(3000,96,artifactDataTime*30000/1000);
    artifactData.t = zeros(3000,1);
    artifactDataIndex = 1;
     
    for i=1:numel(fileList)
        %% load file
        disp(['working on:'])
        disp(fileList(i).name)
        
        cdsTemp=commonDataStructure();
        cdsTemp.file2cds([folderpath,fileList(i).name],inputData.ranBy,inputData.array1,inputData.monkey,inputData.lab,'ignoreJumps',inputData.task,inputData.mapFile,'recoverPreSync');
        
        %% load waveformsSent file if it exists or make one if it does not exist
        underscoreIdx = find(fileList(i).name=='_');
        [~,fname,~] = fileparts(fileList(i).name);
        waveformFilename = strcat(fileList(i).name(1:underscoreIdx(end)),'waveformsSent',fname(underscoreIdx(end):end),'.mat');
        if(exist(waveformFilename)~=0)
            load(waveformFilename)
        end
        
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
%                 error('processStimArtifact:cantFindSync','couldnt find a sync signal')
                noSync = 1;
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
        if(noSync && noSyncIntended)
            
        elseif(noSync)
            artifactDataPre.stimOn = [];
        else
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
        end
        
        %% fix stim times if more than one pulse sent per wave
        if(inputData.moreThanOnePulsePerWave)
            stimOnTemp = [];
            for stimIdx = 1:numel(artifactDataPre.stimOn)
                stimOnTemp = [stimOnTemp; artifactDataPre.stimOn(stimIdx) + 1/inputData.pulseFrequency*30000*(0:1:(inputData.numPulses-1))'];
            end
            artifactDataPre.stimOn = stimOnTemp;
            artifactDataPre.stimOff = artifactDataPre.stimOn + 2;
        end
        %% make waveforms sent file if it does not exist
        if(~exist(waveformFilename)~=0)
            [~,fname,~] = fileparts(fileList(i).name);
            waveformFilename = strcat(fileList(i).name(1:underscoreIdx(end)),'waveformsSent',fname(underscoreIdx(end):end),'.mat');
            waveforms.chanSent = -1*ones(numel(artifactDataPre.stimOn),1);
            waveforms.waveSent = 1*ones(numel(artifactDataPre.stimOn),1);
            waveforms.parameters{1,1} = [];
            save(waveformFilename,'waveforms','-v7.3');
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
                cds.meta.preOffset = preOffset;
                cds.meta.postOffset = postOffset;
                if(~isempty(cdsTemp.kin))
                    cds.kin = cdsTemp.kin;
                end
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
                if(~isempty(cdsTemp.kin))
                    tempKin = cdsTemp.kin;
                    tempKin.t = tempKin.t + cds.meta.duration;
                    cdsKin = cds.kin;
                    cdsKin{end+1:end+size(tempKin,1),:} = tempKin{:,:};
                    cds.kin = cdsKin;
                end
                % port over stimOn and stimOff data
                cds.stimOn = [cds.stimOn;artifactDataPre.stimOn/30000 + cds.meta.duration];
                cds.stimOff = [cds.stimOff;artifactDataPre.stimOff/30000 + cds.meta.duration];
                % port over trial data
                if(~isempty(cdsTemp.trials))
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
                end
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
                if(~isempty(cdsTemp.kin))
                    cds.kin = cdsTemp.kin;
                end
                cds.force = cdsTemp.force;
                cds.lfp = {};
                cds.emg = cdsTemp.emg;
                cds.analog = {};
                if(~noSync)
                    cds.stimOn = artifactDataPre.stimOn/30000;
                    cds.stimOff = artifactDataPre.stimOff/30000;
                end
                cds.triggers = cdsTemp.triggers;
                cds.units = cdsTemp.units;
                cds.aliasList = cdsTemp.aliasList;
                cds.operationLog = cdsTemp.operationLog;
                cds.meta.duration = 0;
            else     
                % port over kin data
                if(~isempty(cdsTemp.kin))
                    tempKin = cdsTemp.kin;
                    tempKin.t = tempKin.t + cds.meta.duration;
                    cdsKin = cds.kin;
                    cdsKin{end+1:end+size(tempKin,1),:} = tempKin{:,:};
                    cds.kin = cdsKin;
                end
                % port over stimOn and stimOff data
                if(~noSync)
                    cds.stimOn = [cds.stimOn;artifactDataPre.stimOn/30000 + cds.meta.duration];
                    cds.stimOff = [cds.stimOff;artifactDataPre.stimOff/30000 + cds.meta.duration];  
                end
                % update meta information
                cds.offsets(end+1,1) = cds.meta.duration;
                cds.meta.dataWindow(2) = cds.meta.dataWindow(2) + cdsTemp.meta.dataWindow(2);
            end
        end
        
        clear cdsTemp
        %% use stim times to filter backwards
        % get template if applicable
        
        if(inputData.templateSubtract && exist(waveformFilename)~=0)
            cdsTempLFPFiltered = cdsTempLFP;
            % need to know stimulated channels
            numStimChan = numel(unique(waveforms.chanSent));
            numStims = zeros(numStimChan,1);
            for stimChanIdx = 1:numStimChan
                stimChans = unique(waveforms.chanSent);
                numStims(stimChanIdx) = sum(waveforms.chanSent == stimChans(stimChanIdx));
            end
            templateAll = zeros(numStimChan,96,templateSize);
            for ch = 2:size(cdsTempLFP,2)
                for stimuli = 1:numel(artifactDataPre.stimOn)+1
                    if(numel(artifactDataPre.stimOn)==0)
                        stimData = cdsTempLFP(:,ch);
                        stimDataTemp = [stimData(:,1);zeros(numZeros,1)];
                        stimDataTemp = fliplr(filter(bFilter,aFilter,fliplr(stimDataTemp')))'; 
                        cdsTempLFPFiltered(:,ch) = stimDataTemp(1:end-numZeros,1); % filter here, override data, move one
                    elseif(stimuli == 1) % all data before first stim
                        stimData = cdsTempLFP(1:artifactDataPre.stimOn(stimuli),ch);
                        stimDataTemp = [stimData(:,1);zeros(numZeros,1)];
                        stimDataTemp = fliplr(filter(bFilter,aFilter,fliplr(stimDataTemp')))'; 
                        cdsTempLFPFiltered(:,ch) = stimDataTemp(1:end-numZeros,1); % filter here, override data, move one
                    elseif(stimuli == numel(artifactDataPre.stimOn) + 1) % all data after last stim artifact
                        stimData = cdsTempLFP(artifactDataPre.stimOn(stimuli-1)+5*30:end,ch);
                        stimDataTemp = [stimData(:,1);zeros(numZeros,1)];
                        stimDataTemp = fliplr(filter(bFilter,aFilter,fliplr(stimDataTemp')))'; 
                        cdsTempLFPFiltered(:,ch) = stimDataTemp(1:end-numZeros,1); % filter here, override data, move one
                    else % data before ith stim up to ith-1 
                        if(artifactDataPre.stimOn(stimuli) - 6*30 > artifactDataPre.stimOn(stimuli-1))  
                            stimData = cdsTempLFP(artifactDataPre.stimOn(stimuli-1):artifactDataPre.stimOn(stimuli),ch);
                            stimChan = find(unique(waveforms.chanSent) == waveforms.chanSent(stimuli-1));
                            stimDataTemp = [stimData(:,1);zeros(numZeros,1)];
                            stimDataTemp = fliplr(filter(bFilter,aFilter,fliplr(stimDataTemp')))';
                            stimDataTemp = stimDataTemp(1:end-numZeros,1);
                            templateAll(stimChan,ch-1,:) = squeeze(templateAll(stimChan,ch-1,:)) + (stimDataTemp(1:templateSize,1)/numStims(stimChan));
                        end
                    end
                    
                end
            end
            
        end
        % get thresholds for each channel based on non stim data
        thresholdAll = zeros(size(cdsTempLFPFiltered,2)-1,1);
        
        for ch = 2:size(cdsTempLFPFiltered,2)
            numPoints = 0;
            for stimuli = 1:numel(artifactDataPre.stimOn)+1
                if(numel(artifactDataPre.stimOn)==0)
                    stimData = cdsTempLFPFiltered(:,ch);
                elseif(stimuli == 1) % all data before first stim
                    stimData = cdsTempLFPFiltered(1:artifactDataPre.stimOn(stimuli),ch);
                elseif(stimuli == numel(artifactDataPre.stimOn) + 1) % all data after last stim artifact
                    stimData = cdsTempLFPFiltered(artifactDataPre.stimOn(stimuli-1)+5*30:end,ch);
                else % data before ith stim up to ith-1 
                    if(artifactDataPre.stimOn(stimuli) - 6*30 > artifactDataPre.stimOn(stimuli-1))  
                        stimData = cdsTempLFPFiltered(artifactDataPre.stimOn(stimuli-1)+6*30:artifactDataPre.stimOn(stimuli),ch);
                    end
                end
%                 numPoints = numPoints + numel(stimData);
%                 stimDataTemp = [stimData(:,1);mean(stimData(end-20:end,1))*ones(numZeros,1)];
%                 stimDataTemp = fliplr(filter(bFilter,aFilter,fliplr(stimDataTemp')))';
%                 
                thresholdAll(ch-1) = thresholdAll(ch-1) + sum(stimData.^2);
            end
            thresholdAll(ch-1) = sqrt(thresholdAll(ch-1)/numPoints);
        end
        
        
        spikeWaves = zeros(10000,lengthWave);
        spikeTimes = zeros(10000,1);
        spikeChan = zeros(10000,1);    
        spikeNum = 1;
        disp('filtering and thresholding')
        channelsStimulated = unique(waveforms.chanSent);
        for channelSent = 1:numel(channelsStimulated)
            disp(num2str(channelSent))
            for waveSent = 1:numel(waveforms.parameters)
                disp(num2str(waveSent))
                for ch = (1:size(cdsTempLFPFiltered,2)-1)      
%                     ch=2;
                    %% build stim data matrix for this channel
                    stimData = [];
                    stimDataSizes = [];
                    for stimIdx = 1:numel(artifactDataPre.stimOn)+1
                        if((stimIdx == 1 && channelSent == 1 && waveSent == 1) || ... 
                                (stimIdx~=1 && waveforms.chanSent(stimIdx-1) == channelsStimulated(channelSent) && waveforms.waveSent(stimIdx-1) == waveSent))
                            if(numel(artifactDataPre.stimOn) == 0) % no stimulations
                                stimDataTemp = cdsTempLFPFiltered(1:inputData.artifactDataTime*30,ch+1);
                            elseif(stimIdx == 1) % all data before first stim -- do nothing since no stimulation
%                                 stimDataTemp = cdsTempLFP(1:artifactDataPre.stimOn(stimIdx),ch+1);
                            elseif(stimIdx == numel(artifactDataPre.stimOn) + 1) % all data after last stim
                                stimDataTemp = cdsTempLFPFiltered(artifactDataPre.stimOn(stimIdx-1):artifactDataPre.stimOn(stimIdx-1)+inputData.artifactDataTime*30,ch+1);
                %                 % perform pca step
                %                 [coeff,score] = pca(stimData);
                %                 coeff(:,1:4) = 0;
                %                 stimData = coeff*score';
                %                 stimData = stimData';
            %                     if(~inputData.templateSubtract)
            %                         stimData(1:30*1,:) = 0; % blank the first millisecond
            %                     end
                            else % data before ith stim up to ith-1 
                                stimDataTemp = cdsTempLFPFiltered(artifactDataPre.stimOn(stimIdx-1):artifactDataPre.stimOn(stimIdx-1)+inputData.artifactDataTime*30,ch+1);
                                % perform pca step
                %                 [coeff,score] = pca(stimData);
                %                 coeff(:,1:4) = 0;
                %                 stimData = coeff*score';
                %                 stimData = stimData';
            %                     if(~inputData.templateSubtract)
            %                         stimData(1:30*1,:) = 0; % blank the first millisecond
            %                     end
                            end

                            % join stimDataTemp with stimData
                            stimDataTemp = stimDataTemp';
                            stimDataSizes(end+1,1) = size(stimDataTemp,2);
                            if(stimIdx == 1 || isempty(stimData))
                                stimData = stimDataTemp;
                            else
                                stimData(end+1,:) = stimDataTemp;
                            end
%                             if(stimIdx == 1 || isempty(stimData))
%                                 stimData = stimDataTemp;
%                             elseif(size(stimDataTemp,2) > size(stimData,2))
%                                 numPad = size(stimDataTemp,2) - size(stimData,2);
%                                 for numPadIdx = 1:numPad
%                                     stimData(:,end+1) = mean(stimData(:,end-10:end),2);
%                                 end
%                                 stimData(end+1,:) = stimDataTemp;
%                             elseif(size(stimDataTemp,2) < size(stimData,2))
%                                 numPad = size(stimData,2) - size(stimDataTemp,2);
%                                 stimDataTemp(1,end+1:end+numPad) = mean(stimDataTemp(:,end-10:end));
%                                 stimData(end+1,:) = stimDataTemp;
%                             else
%                                 stimData(end+1,:) = stimDataTemp;   
%                             end
                        end
                    end


                    %% filter backwards for all stimulations on this channel -- already did this step above
%                     stimDataPad = stimData;
%                     stimDataPad(:,end+1:end+numZeros) = repmat(mean(stimDataPad,2),1,numZeros);
%                     stimDataFilt = fliplr(filter(bFilter,aFilter,fliplr(stimDataPad)')');
%                     stimDataFilt = stimDataFilt(:,1:end-numZeros);
%                     clear stimDataPad
%                     clear stimData
                    
                    %% spline interpolate for all stimulations on this channel
                    xFilt = 1:1:size(stimDataFilt,2);
                    xQuery = 1:1/inputData.interpolateValue:size(stimDataFilt,2);
                    stimDataSplined = zeros(size(stimDataFilt,1),inputData.artifactDataTime*30*inputData.interpolateValue + size(stimDataFilt,2) - inputData.artifactDataTime*30);
                    for a = 1:size(stimDataFilt,1)
                        splineTemp = spline(xFilt(1:inputData.artifactDataTime*30),stimDataFilt(a,1:inputData.artifactDataTime*30),xQuery(1:inputData.artifactDataTime*30*inputData.interpolateValue));
                        stimDataSplined(a,:) = [splineTemp,stimDataFilt(a,inputData.artifactDataTime*30+1:end)];
                    end
                    clear stimDataFilt
                    
                    %% find template
                    template = mean(stimDataSplined);

                    %% shift in time
                    resGain = inputData.interpolateValue;
                    resBound = 5;
                    stimDataShifted = stimDataSplined;
                    clear stimDataSplined
                    sse = [];
                    shift = zeros(size(stimDataShifted,1),1);
                    for w = 1:size(stimDataShifted,1)
                        % test every possible shift value, pick best
                        for shiftValue = (-resGain*resBound+1):1:(resGain*resBound)
                            testWave = stimDataShifted(w,1:inputData.artifactDataTime*30*inputData.interpolateValue);
                            testWave = circshift(testWave,shiftValue);
                            if(shiftValue < 0)
                                testWave(end+shiftValue:end) = 0;
                            else
                                testWave(1:shiftValue) = 0;
                            end
                            sse(shiftValue+resGain*resBound) = sum((testWave-template(1:numel(testWave))).^2);
                        end
                        sse(sse==0) = max(sse);
                        [~,minIdx] = min(sse);
                        shift(w) = minIdx-resGain*resBound;
                        stimDataShifted(w,:) = circshift(stimDataShifted(w,:),shift(w));
                        if(shift(w) < 0)
                            stimDataShifted(end+shift(w):end) = 0;
                        else
                            stimDataShifted(1:shift(w)) = 0;
                        end
                    end
                    %% scale amplitude
                    stimDataScaled = stimDataShifted;
        %             clear stimDataShifted
                    sse = [];
                    scaleValues = (0.8:0.01:1.2)';
                    for w = 1:size(stimDataScaled,1)
                        % test every possible shift value, pick best
                        testWave = scaleValues*stimDataScaled(w,1:inputData.artifactDataTime*30*inputData.interpolateValue);
                        sse = sum((testWave-repmat(template(1,1:size(testWave,2)),size(testWave,1),1)).^2');
                        [~,minIdx] = min(sse);
                        scale(w) = scaleValues(minIdx);
                        stimDataScaled(w,:) = scaleValues(minIdx)*stimDataScaled(w,:);
                    end
                    %% find template
                    template = mean(stimDataScaled);

                    %% subtract template
                    stimDataAll = stimDataScaled - template;

                    %% unscale
                    try
                        stimDataAll = stimDataAll./scale';
                    catch
                        disp('Warning: something minor happened');
                    end
                    %% unshift
                    for a = 1:size(stimDataAll,1)
                        stimDataAll(w,:) = circshift(stimDataAll(w,:),-1*shift(w));
                        if(-1*shift(w) < 0)
                            stimDataAll(end+shift(w):end) = 0;
                        else
                            stimDataAll(1:shift(w)) = 0;
                        end
                    end

                    %% downsample
                    stimDataAll = [downsample(stimDataAll(:,1:inputData.artifactDataTime*30*inputData.interpolateValue)',inputData.interpolateValue)',stimDataAll(:,inputData.artifactDataTime*30*inputData.interpolateValue+1:end)];

                    %% blank
                    stimDataAll(:,1:inputData.templateBlankPeriod) = 0;

                    %% for each stim idx
                    stimDataIdx = 0;
                    for stimIdx = 1:numel(artifactDataPre.stimOn)+1
                        if((stimIdx == 1 && channelSent == 1 && waveSent == 1) || ... 
                                (stimIdx~=1 && waveforms.chanSent(stimIdx-1) == channelsStimulated(channelSent) && waveforms.waveSent(stimIdx-1) == waveSent))
                            stimDataIdx = stimDataIdx + 1;
                            
                            %% add back data removed previously
                            if(numel(artifactDataPre.stimOn) == 0) % no stimulations
                                stimData = cdsTempLFPFiltered(:,ch+1);
                            elseif(stimIdx == 1) % all data before first stim -- grab all data
                                stimData = cdsTempLFPFiltered(1:artifactDataPre.stimOn(stimIdx),ch+1);
                            elseif(stimIdx == numel(artifactDataPre.stimOn) + 1) % all data after last stim
                                stimData = [stimDataAll(stimDataIdx,:)';cdsTempLFPFiltered(artifactDataPre.stimOn(stimIdx-1)+inputData.artifactDataTime*30+1:end,ch+1)];
                                cdsTempLFPFiltered(artifactDataPre.stimOn(stimIdx-1):artifactDataPre.stimOn+inputData.artifactDataTime*30,ch+1) = stimDataAll(stimDataIdx,:)';
                            else % data before ith stim up to ith-1 
                                stimData = [stimDataAll(stimDataIdx,:)';cdsTempLFPFiltered(artifactDataPre.stimOn(stimIdx-1)+inputData.artifactDataTime*30+1:artifactDataPre.stimOn(stimIdx),ch+1)];
                                cdsTempLFPFiltered(artifactDataPre.stimOn(stimIdx-1):artifactDataPre.stimOn(stimIdx-1)+inputData.artifactDataTime*30,ch+1) = stimDataAll(stimDataIdx,:)';
                            end
%                             %% remove excess data from joining
%                             stimData = stimDataAll(stimDataIdx,1:stimDataSizes(stimDataIdx));
%                             stimData = stimData';    
                           %% get threshold and threshold crossings
            %                 threshold = abs(rms(stimData(max(1,numel(stimData(:,ch))-10):end,ch))*thresholdMult);
                            threshold = inputData.thresholdMult*thresholdAll(ch);
                            thresholdCrossings = find(stimData(:,1)>threshold);

                            % remove potential artifacts and too close to beginning
                            crossingsMask = ones(numel(thresholdCrossings),1);
                            for cross = 1:numel(thresholdCrossings)
                                if(stimData(thresholdCrossings(cross),1) > maxAmplitude || ...
                                        ~(thresholdCrossings(cross)+postOffset <= numel(stimData(:,1)) && thresholdCrossings(cross)-preOffset > 0))
                                    crossingsMask(cross) = 0;
                                end
                            end
                            thresholdCrossings = thresholdCrossings(crossingsMask(:) == 1);


                            % remove chains -- find best spot
                            idx = 2;
                            chain = [1];
                            crossingsKeep = [];
                            while idx <= numel(thresholdCrossings)
                                if(thresholdCrossings(idx) == thresholdCrossings(idx-1)+1) % store in chain
                                    chain = [chain;idx];
                                elseif(~isempty(chain)) % broke a chain, store minidx, update idx, empty chain
                                    [~,maxIdx] = max(stimData(thresholdCrossings(chain)));
            %                         [~,maxIdx] = max(chain);
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
                                spikeWaves(spikeNum,:) = stimData(thresholdCrossings(cross)-preOffset:thresholdCrossings(cross)+postOffset,1);
                                spikeNum = spikeNum + 1;

                                if(spikeNum > 0.67*numel(spikeTimes))
                                    spikeTimes = [spikeTimes;zeros(1000,1)];
                                    spikeWaves = [spikeWaves; zeros(1000,lengthWave)];
                                    spikeChan = [spikeChan; zeros(1000,1)];
                                end
                            end
                        end
                    end % stim idx for
               end % end channel for
            end
        end
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
                artifactData.artifactProcessed(artifactDataIndex,:,:) = cdsTempLFPFiltered(artifactDataPre.stimOn(art):artifactDataPre.stimOn(art)+floor(artifactDataTime*30000/1000)-1,2:end)';
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
    cds.units = [];
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

