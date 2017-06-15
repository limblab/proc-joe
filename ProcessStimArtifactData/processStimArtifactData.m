function [outputFigures, outputData ] = processStimArtifact(folderpath, inputData )
    %script to load stimulation files and generate perievent plots of 30khz
    %data. Formatted to work with runDataProcessing
    outputFigures = [];
    outputData=[];
    %get list of all files in the folder:
    if ~strcmp(folderpath(end),filesep)
        folderpath=[folderpath,filesep];
    end
    cd(folderpath);
    fileList=dir('*001.nev');

    mapData=loadMapFile(inputData.mapFile(8:end));
        for i=1:size(mapData,1)
            mapData.label{i}=[inputData.array1(6:end),mapData.label{i}];
        end
    %establish trackers for information that is common to all files:
    eList={};
    posList=[];
    chList=[];
    for i=1:numel(fileList)
        %% load file
        disp(['working on:'])
        disp(fileList(i).name)
        
        cds=commonDataStructure();
        cds.file2cds([folderpath,fileList(i).name],inputData.ranBy,inputData.array1,inputData.monkey,inputData.lab,'ignoreJumps',inputData.task,inputData.mapFile,'recoverPreSync');
        
        %% find sync signal in analog data
        useSync=true;
        aIdx=[];
        syncIdx=[];
        if isempty(inputData.useSyncLabel)
        %look for sync under the label sync
            for j=1:numel(cds.analog)
                syncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'sync'));
                if ~isempty(syncIdx)
                    aIdx=j;
                    syncName='sync';
                end
            end
            %if it wasn't called sync, try for matt's 'StimTrig' label:
            if isempty(syncIdx)
                for j=1:numel(cds.analog)
                    syncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'StimTrig'));
                    if ~isempty(syncIdx)
                        aIdx=j;
                        syncName='StimTrig';
                    end
                end
            end
            %if we didn't find a sync channel, just look for ainp16
            if isempty(syncIdx)
                for j=1:numel(cds.analog)
                    syncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'ainp16'));
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
                for j=1:numel(cds.analog)
                    syncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'sync'));
                    if ~isempty(syncIdx)
                        aIdx=j;
                        syncName='sync';
                    end
                end
            else
                for j=1:numel(cds.analog)
                    syncIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,'ainp16'));
                    if ~isempty(syncIdx)
                        useSync=false;
                        aIdx=j;
                        syncName='ainp16';
                    end
                end
            end
        end
        %% use sync to get stim times:
        artifactData(i).stimOn=find(diff(cds.analog{aIdx}.(syncName)-mean(cds.analog{aIdx}.(syncName))>3)>.5);
        stimOff=find(diff(cds.analog{aIdx}.(syncName)-mean(cds.analog{aIdx}.(syncName))<-3)>.5);
        artifactData(i).stimOn = artifactData(i).stimOn;
        artifactData(i).stimOff=nan(size(artifactData(i).stimOn));
        for j=1:numel(artifactData(i).stimOn)
            if j<numel(artifactData(i).stimOn)
                next=artifactData(i).stimOn(j+1);
            else
                next=numel(cds.analog{aIdx}.(syncName));
            end
            offIdx=stimOff(find((stimOff>artifactData(i).stimOn(j)& stimOff<next),1,'first'));
            if ~isempty(offIdx)
                artifactData(i).stimOff(j)=offIdx;
            end
        end

        stimWindows=[artifactData(i).stimOn,artifactData(i).stimOn+inputData.windowSize-1];    
        %% find monitor data for cerestim monitoring ports:
        monitorName='monitor';
        for j=1:numel(cds.analog)
            monitorIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,[monitorName,'1']));
            if ~isempty(monitorIdx)
                mIdx=j;
            end
        end

        %% put spikes into structure:
        idxStart=strfind(fileList(i).name,'chan')+4;
        if numel(idxStart)>1
            %we have 2 instances of chan in the name, try again looking
            %for _chan
            idxStart=strfind(fileList(i).name,'_chan')+5;
        end
        idxEnd=strfind(fileList(i).name,'stim')-1;
        if numel(idxEnd)>1
            %we have 2 instances of end, look for stim_ instead and
            %hope that clears it up:
            idxEnd=strfind(fileList(i).name,'stim_')-1;
        end
        if isempty(idxStart)
            %hopefully this is a file that has the ch#stim format:
            for j=1:numel(idxEnd)
                idxStart=strfind(fileList(i).name(1:idxEnd(j)),'ch');
                if isempty(idxStart)
                    continue
                end
                idxStart=idxStart(end)+2;
                if idxEnd(j)-idxStart<4
                    idxEnd=idxEnd(end);
                end
            end
        end
        artifactData(i).stimChannel=str2num(fileList(i).name(idxStart:idxEnd));

        if isempty(cds.lfp)
            numChans=size(cds.analog{1,1},2);
        else
            numChans=size(cds.lfp,2);
        end
        
        
        artifactMat=[];
        monitorCh1Mat=[];
        monitorCh2Mat=[];
        for j=2:numChans
            for k=1:numel(artifactData(i).stimOn)
                %disp(['j=',num2str(j),' k=',num2str(k)])
                if isempty(cds.lfp)
                    artifactMat(j-1,k,1:inputData.windowSize+inputData.presample)=reshape(cds.analog{1,1}{stimWindows(k,1)-inputData.presample:stimWindows(k,2),j},[1,1,inputData.windowSize+inputData.presample]);
                    chanNum=str2num(cds.analog{1,1}.Properties.VariableNames{j}(5:end));
                    mapIdx=find(mapData.chan==chanNum);
                    electrodeList{j-1}=mapData.label(mapIdx);
                else
                    artifactMat(j-1,k,1:inputData.windowSize+inputData.presample)=reshape(cds.lfp{stimWindows(k,1)-inputData.presample:stimWindows(k,2),j},[1,1,inputData.windowSize+inputData.presample]);
                    electrodeList{j-1}=cds.lfp.Properties.VariableNames{j};
                end
                %if we don't have a position for this electrode, find it:
                if isempty(find(strcmp(eList,electrodeList{j-1}),1))
                    %find the electrode in the units structure and add it
                    %to eList and posList:
                    unitIdx=find(strcmp(mapData.label,electrodeList{j-1}));
                    if ~isempty(unitIdx)
                        eList(end+1)=electrodeList(j-1);
                        chList(end+1)=(uint8(mapData.bank{unitIdx})-65)*32+mapData.pin(unitIdx);
                        posList(end+1,:)=[mapData.row(unitIdx),mapData.col(unitIdx)];
                    end
                end
                if exist('mIdx','var')
                    monitorCh1Mat(j-1,k,1:inputData.windowSize+inputData.presample)=reshape(cds.ananlog{mIdx}.([monitorName,'1'])(stimWindows(k,1)-inputData.presample:stimWindows(k,2)),[1,1,inputData.windowSize+inputData.presample]);
                    monitorCh2Mat(j-1,k,1:inputData.windowSize+inputData.presample)=reshape(cds.ananlog{mIdx}.([monitorName,'2'])(stimWindows(k,1)-inputData.presample:stimWindows(k,2)),[1,1,inputData.windowSize+inputData.presample]);
                end
            end
        end
        artifactData(i).artifact=artifactMat;
        artifactData(i).electrodeNames=electrodeList;
        artifactData(i).monitor1=monitorCh1Mat;
        artifactData(i).monitor2=monitorCh2Mat;
        %now get the monitor data showing the driving voltage the
        %stimulator used for the pulses
        for j=2:numChans
            for k=1:numel(artifactData(i).stimOn)
            end
        end


        %clear cds
        clear cds
    end
    outputData.artifactData=artifactData;
    outputData.eList=eList;
    outputData.posList=posList;
    outputData.chList=chList;
 
end
