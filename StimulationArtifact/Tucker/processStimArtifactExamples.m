function [ outputFigures,outputData ] = processStimArtifactExamples(folderpath, inputData )
    %script to load stimulation files and generate perievent plots of 30khz
    %data. Formatted to work with runDataProcessing
    
    outputFigures=[];
    outputData=[];
    %get list of all files in the folder:
    if ~strcmp(folderpath(end),filesep)
        folderpath=[folderpath,filesep];
    end
    cd(folderpath);
    fileList=dir('*.nev');

    %get the mapfile info:
%     mapFile=inputData.mapFile(8:end);
%     mapData=readtable(mapFile,'FileType','text','HeaderLines',13,'Delimiter','tab');
%     mapData.Properties.VariableNames{1}='col';%fix the column header
%     %fix the labels:
%     for i=1:size(mapData,1)
%         mapData.label{i}=[inputData.array1(6:end),mapData.label{i}];
%     end
    
    mapData=loadMapFile(inputData.mapFile(8:end));
        for i=1:size(mapData,1)
            mapData.label{i}=[inputData.array1(6:end),mapData.label{i}];
        end
    %establish trackers for information that is common to all files:
    eList={};
    posList=[];
    chList=[];
    if RDPIsAlreadyDone('artifactData',folderpath)
        warning('processStimArtifact:foundExistingData','loading data from previous processing. This will have the PREVIOUS settings for time window, presample etc')
        artifactData =RDPLoadExisting('artifactData',folderpath);
        eList =RDPLoadExisting('eList',folderpath);
        posList =RDPLoadExisting('posList',folderpath);
        chList =RDPLoadExisting('chList',folderpath);
    else
        for i=1:numel(fileList)
            %% load file
            disp(['working on:'])
            disp(fileList(i).name)
            cds=commonDataStructure();
            cds.file2cds([folderpath,fileList(i).name],inputData.ranBy,inputData.array1,inputData.monkey,inputData.lab,'ignoreJumps',inputData.task,inputData.mapFile,inputData.useBlock);
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
            stimOn=find(diff(cds.analog{aIdx}.(syncName)-mean(cds.analog{aIdx}.(syncName))>3)>.5);
%             if useSync
%                 stimOn=find(diff(cds.analog{aIdx}.sync-mean(cds.analog{aIdx}.sync)>3)>.5);
%             else
%                 stimOn=find(diff(cds.analog{aIdx}.ainp16-mean(cds.analog{aIdx}.ainp16)>3)>.5);
%             end
%             if numel(stimOn>100)
%                 %loop through till we find a likely channel:
%                 for i=2:size(cds.lfp,2)
%                 end
%             end
            stimWindows=[stimOn,stimOn+inputData.windowSize-1];    
            %put spikes into structure:
            idxStart=strfind(fileList(i).name,'chan')+4;
            idxEnd=strfind(fileList(i).name,'stim')-1;
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

%             labelList={cds.units.label};
%             arrayList={cds.units.array};
%             for j=1:numel(labelList)
%                 labelList{j}=[arrayList{j},labelList{j}];
%             end
            for j=2:size(cds.lfp,2)
                for k=1:numel(stimOn)
                    electrodeList{j-1}=cds.lfp.Properties.VariableNames{j};
                    artifactMat(j-1,k,1:inputData.windowSize+inputData.presample)=reshape(cds.lfp{stimWindows(k,1)-inputData.presample:stimWindows(k,2),j},[1,1,inputData.windowSize+inputData.presample]);
                    %if we don't have a position for this electrode, find it:
                    if isempty(find(strcmp(eList,cds.lfp.Properties.VariableNames{j}),1))
                        %find the electrode in the units structure and add it
                        %to eList and posList:
                        unitIdx=find(strcmp(mapData.label,cds.lfp.Properties.VariableNames{j}));
                        if ~isempty(unitIdx)
                            eList(end+1)=cds.lfp.Properties.VariableNames(j);
                            chList(end+1)=(uint8(mapData.bank{unitIdx})-65)*32+mapData.pin(unitIdx);
                            posList(end+1,:)=[mapData.row(unitIdx),mapData.col(unitIdx)];
                        end
                    end
                end
            end
            artifactData(i).artifact=artifactMat;
            artifactData(i).electrodeNames=electrodeList;
            %clear cds
            clear cds
        end
        outputData.artifactData=artifactData;
        outputData.eList=eList;
        outputData.posList=posList;
        outputData.chList=chList;
    end
    
    %plot the artifacts for each stim channel:
    outputFigures(end+1)=figure;
    set(outputFigures(end),'Name','allStimChans')
    numPlotPixels=1200;
    set(outputFigures(end),'Position',[100 100 numPlotPixels numPlotPixels]);
    paperSize=0.2+numPlotPixels/get(outputFigures(end),'ScreenPixelsPerInch');
    set(outputFigures(end),'PaperSize',[paperSize,paperSize]);
    timevec=([1:size(artifactData(1).artifact,3)]-inputData.presample-1)/30;
    for i=2:numel(artifactData)
        for j=1:numel(artifactData(i).electrodeNames)
        %    disp(['   plotting response for electrode:',num2str(j)])
            %find the index in our full list of electrodes, that
            %corresponds to the channel j in artifactData(i)
            posIdx=find(strcmp(eList,artifactData(i).electrodeNames{j}));
            if chList(posIdx)==artifactData(i).stimChannel
                %set color for main artifact plot
                colorCode='r';
                figure(outputFigures(end));
                subplot(numel(artifactData)-1,1,i-1);
                %plot(timevec, squeeze(artifactData(i).artifact(j,:,:))'-squeeze(repmat(artifactData(i).artifact(j,:,1),[1,1,size(artifactData(i).artifact,3)]))',colorCode)
                plot(timevec, squeeze(artifactData(i).artifact(j,:,:))',colorCode)
                %add line @ 1ms:
                hold on
                plot([1,1]*(inputData.presample+30)/30,[8000,-8000],'g')
                %add line @ 1.5ms:
                plot([1,1]*(inputData.presample+45)/30,[8000,-8000],'b')
                axis tight
            end
        end
    end
    set(gcf,'NextPlot','add');
    axes;
    h = title(['Artifact on stim Channels ']);
    set(gca,'Visible','off');
    set(h,'Visible','on'); 
    
    outputFigures(end+1)=figure;
    set(outputFigures(end),'Name','randomChans')
    numPlotPixels=1200;
    set(outputFigures(end),'Position',[100 100 numPlotPixels numPlotPixels]);
    paperSize=0.2+numPlotPixels/get(outputFigures(end),'ScreenPixelsPerInch');
    set(outputFigures(end),'PaperSize',[paperSize,paperSize]);
    chanList=[4 73];
    for i=2:numel(artifactData)
        for j=1:numel(chanList)
            %set color for main artifact plot
            colorCode='k';
            figure(outputFigures(end));
            subplot(numel(artifactData)-1,numel(chanList),(i-2)*numel(chanList)+j);
            %figure
            plot(timevec,squeeze(artifactData(i).artifact(chanList(j),:,:))'-squeeze(repmat(artifactData(i).artifact(chanList(j),:,1),[1,1,size(artifactData(i).artifact,3)]))',colorCode)
            %add line @ 1ms:
            hold on
            plot([1,1]*(inputData.presample+30)/30,[8000,-8000],'g')
            %add line @ 1.5ms:
            plot([1,1]*(inputData.presample+45)/30,[8000,-8000],'b')
            axis tight
            ylim([-inputData.plotRange,inputData.plotRange]*1000);
        end
    end
    set(gcf,'NextPlot','add');
    axes;
    h = title(['Artifact in response to stim on Channels: ',num2str(chanList)]);
    set(gca,'Visible','off');
    set(h,'Visible','on'); 
end

