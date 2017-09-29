function [artifact] = extractArtifactData(folderpath, inputData)

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
    thresholdMult = inputData.thresholdMult;
    % variables to store spike information   
%     preOffset = 27;
%     postOffset = 20;
    preOffset = inputData.preOffset;
    postOffset = inputData.postOffset;
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

        cdsTemp.file2cds([folderpath,fileList(i).name],inputData.ranBy,inputData.array1,inputData.monkey,inputData.lab,inputData.task,inputData.mapFile);
        
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
            artifactDataPre.stimOn = [];
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
        if(~exist(waveformFilename)~=0 && ~(noSync || noSyncIntended))
            [~,fname,~] = fileparts(fileList(i).name);
            waveformFilename = strcat(fileList(i).name(1:underscoreIdx(end)),'waveformsSent',fname(underscoreIdx(end):end),'.mat');
            waveforms.chanSent = -1*ones(numel(artifactDataPre.stimOn),1);
            waveforms.waveSent = 1*ones(numel(artifactDataPre.stimOn),1);
            waveforms.parameters{1,1} = [];
            save(waveformFilename,'waveforms','-v7.3');
        end
        %% extract data from cdsTemp so that the rest of the code moves quicker - table operations are slow af
        cdsTempLFP = cdsTemp.lfp{:,:};
        if(cdsTempLFP(end,1) > cdsTemp.meta.duration-0.01) % apparently need to remove data?
            timeRemove = cdsTempLFP(end,1) - cdsTemp.meta.duration + 0.2;
            pointsRemove = timeRemove*30000;
            cdsTempLFP = cdsTempLFP(1:end-pointsRemove,:);
        end
        cdsTempDuration = cdsTemp.meta.duration;
        
        artifact = zeros(numel(artifactDataPre.stimOn)-1,32,1000);
        for st = 1:numel(artifactDataPre.stimOn)-1
            artifact(st,:,:) = squeeze(cdsTempLFP(artifactDataPre.stimOn(st)-99:artifactDataPre.stimOn(st)+900,2:33)');
        end
    end


end