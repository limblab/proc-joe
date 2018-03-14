function [ outputFigures,outputData ] = processStimArtifact_filter(folderpath, inputData )
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
    
    mapData=loadMapFile(inputData.mapFile(8:end));
        for i=1:size(mapData,1)
            mapData.label{i}=[inputData.array1(6:end),mapData.label{i}];
        end
    %establish trackers for information that is common to all files:
    eList={};
    posList=[];
    chList=[];
    interleavedTrials = 0;
    
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
%         artifactData(i).stimOn = artifactData(i).stimOn(1);
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
            monitorIdx=find(strcmp(cds.analog{j}.Properties.VariableNames,[monitorName,'1']) | strcmp(cds.analog{j}.Properties.VariableNames,'Cerestim_module'));
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
            numChans=size(cds.analog{1,1},2)-1;
        else
            numChans=size(cds.lfp,2)-1;
        end
        
        % check for duke board and replace lfp data with the analog data
        if(inputData.dukeBoardChannel < 1 && (~isempty(strfind(fileList(i).name,'dukeBoard')) || ~isempty(strfind(fileList(i).name,'chic201802')))) % search for dukeBoard or chic201802 in the name
            inputData.dukeBoardChannel = artifactData(i).stimChannel;
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
                elseif(j-1 == inputData.dukeBoardChannel)
                    artifactMat(j-1,k,1:inputData.windowSize+inputData.presample)=reshape(cds.analog{1,1}.(inputData.dukeBoardLabel)(stimWindows(k,1)-inputData.presample:stimWindows(k,2)),[1,1,inputData.windowSize+inputData.presample]);
                    electrodeList{j-1}=cds.lfp.Properties.VariableNames{j};
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
                    monitorCh1Mat(j-1,k,1:inputData.windowSize+inputData.presample)=reshape(cds.analog{mIdx}.Cerestim_module(stimWindows(k,1)-inputData.presample:stimWindows(k,2)),[1,1,inputData.windowSize+inputData.presample]);
%                     monitorCh1Mat(j-1,k,1:inputData.windowSize+inputData.presample)=reshape(cds.ananlog{mIdx}.([monitorName,'1'])(stimWindows(k,1)-inputData.presample:stimWindows(k,2)),[1,1,inputData.windowSize+inputData.presample]);
%                     monitorCh2Mat(j-1,k,1:inputData.windowSize+inputData.presample)=reshape(cds.ananlog{mIdx}.([monitorName,'2'])(stimWindows(k,1)-inputData.presample:stimWindows(k,2)),[1,1,inputData.windowSize+inputData.presample]);
                end
            end
        end
        
        artifactData(i).artifact=artifactMat;
        artifactData(i).electrodeNames=electrodeList;
        artifactData(i).monitor1=monitorCh1Mat;
        artifactData(i).monitor2=monitorCh2Mat;
        
        % filter the data in artifactData(i).artifact and store as
        % artifactFiltered
        artifactData(i).artifactFiltered = artifactData(i).artifact;
        for chanIdx = 1:size(artifactData(i).artifact,1)
            artifactData(i).artifactFiltered(chanIdx,:,:) = acausalFilter(squeeze(artifactData(i).artifact(chanIdx,:,:))')';
        end
        
        % look for a 'waveformsSent' file in the directory 
        waveformsSentFileList = dir('*waveformsSent*');
        waveformsIdx = [];
        for w = 1:numel(waveformsSentFileList)
            % remove _waveformsSent from file name, compare to file name
            % currently on
            nameIdx = strfind(waveformsSentFileList(w).name,'waveformsSent');
            waveformsSentFileName = strcat(waveformsSentFileList(w).name(1:nameIdx-2),waveformsSentFileList(w).name(nameIdx+13:end));
            if(strcmpi(waveformsSentFileName(1:end-4),fileList(i).name(1:end-4)))
                waveformsIdx = w;
                interleavedTrials = 1;
            end
        end
        
        artifactData(i).waveNum = zeros(size(artifactData(i).stimOn));
        artifactData(i).chanSent = zeros(size(artifactData(i).stimOn));
        if(~isempty(waveformsIdx))
            load(waveformsSentFileList(waveformsIdx).name); % load in file
            artifactData(i).waveNum = waveforms.waveSent(1:numel(artifactData(i).stimOn));
            artifactData(i).chanSent = waveforms.chanSent(1:numel(artifactData(i).stimOn));
            artifactData(i).parameters = waveforms.parameters;
        else % no waveformsSent file, assume cathodal then anodal alternating pulses
            artifactData(i).waveNum(1:2:numel(artifactData(i).stimOn)) = 1;
            artifactData(i).waveNum(2:2:numel(artifactData(i).stimOn)) = 2;
            artifactData(i).chanSent = ones(size(artifactData(i).stimOn))*artifactData(i).stimChannel;
            artifactData(i).parameters = [];
        end
        
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
    

    %% plot the artifacts for each stim channel/waveform:
    if(interleavedTrials) % only plot data from artifactData(1), but make a 
        % plot for each waveform/stim channel combination
        numChans = numel(unique(artifactData(1).chanSent));
        stimChanList = unique(artifactData(1).chanSent);
        numWaves = numel(unique(artifactData(1).waveNum));
        for nc = 1:numChans
            for nw = 1:numWaves
                outputFigures(end+1) = figure;
                set(outputFigures(end),'Name',['StimOnCH-',num2str(stimChanList(nc)),'Wave-',num2str(nw),'_AllChannels'])
                numPlotPixels=1200;
                set(outputFigures(end),'Position',[100 100 numPlotPixels numPlotPixels]);
                paperSize=0.2+numPlotPixels/get(outputFigures(end),'ScreenPixelsPerInch');
                set(outputFigures(end),'PaperSize',[paperSize,paperSize]);
                for j = 1:numel(artifactData(1).electrodeNames)
                    %find the index in our full list of electrodes, that
                    %corresponds to the channel j in artifactData(i)
                    posIdx=find(strcmp(eList,artifactData(1).electrodeNames{j}));
                    eRow=posList(posIdx,1);
                    eCol=posList(posIdx,2);
                    h=subplot(10,10,10*(eRow-1)+eCol);
                    hold on
                    if chList(posIdx)==stimChanList(nc)
                        %put a purple box around the stim channel:
                        stimBoxX=[.95*(-inputData.presample/30),.95*(size(artifactData(1).artifact,3)-inputData.presample)/30,.95*(size(artifactData(1).artifact,3)-inputData.presample)/30,...
                            .95*(-inputData.presample/30),.95*(-inputData.presample/30)];
                        stimBoxY=[-inputData.plotRange*.95*1000,-inputData.plotRange*.95*1000,inputData.plotRange*.95*1000,inputData.plotRange*.95*1000,-inputData.plotRange*.95*1000];
                        plot(stimBoxX,stimBoxY,'mp-','linewidth',3)
                    end
                    if find(chList(posIdx)==inputData.badChList,1,'first')
                       %put a red X through known bad channels:
                       badMarkX=[2,.95*size(artifactData(1).artifact,3)];
                       badMarkY=[-inputData.plotRange*.95*1000,inputData.plotRange*.95*1000];
                       plot(badMarkX,badMarkY,'rp-','lineWidth',3)
                       badMarkX=[2,.95*size(artifactData(1).artifact,3)];
                       badMarkY=[inputData.plotRange*.95*1000,-inputData.plotRange*.95*1000];
                       plot(badMarkX,badMarkY,'rp-','lineWidth',3)
                    end
                    inputData.filtered = 0;
                    inputData.plotAllChannels = 1;
                    plotArtifacts_singleParameter(inputData,artifactData,j,stimChanList(nc),nw);
                end
            end
        end
        
        % make non-filtered plot of stimulated channel
        for nc = 1:numChans
            for nw = 1:numWaves
                outputFigures(end+1) = figure;
                set(outputFigures(end),'Name',['StimOnCH-',num2str(stimChanList(nc)),'Wave-',num2str(nw),'_Raw'])
                for j = 1:numel(artifactData(1).electrodeNames)
                    %find the index in our full list of electrodes, that
                    %corresponds to the channel j in artifactData(i)
                    posIdx=find(strcmp(eList,artifactData(1).electrodeNames{j}));
                    hold on
                    if chList(posIdx)==stimChanList(nc)
                        % make plot
                        inputData.filtered = 0;
                        inputData.plotAllChannels = 0;  
                        plotArtifacts_singleParameter(inputData,artifactData,j,stimChanList(nc),nw);
                    end
                end
            end
        end
        
        % make filtered plot of stimulated channel
        for nc = 1:numChans
            for nw = 1:numWaves
                outputFigures(end+1) = figure;
                set(outputFigures(end),'Name',['StimOnCH-',num2str(stimChanList(nc)),'Wave-',num2str(nw),'_Filtered'])
                for j = 1:numel(artifactData(1).electrodeNames)
                    %find the index in our full list of electrodes, that
                    %corresponds to the channel j in artifactData(i)
                    posIdx=find(strcmp(eList,artifactData(1).electrodeNames{j}));
                    hold on
                    if chList(posIdx)==stimChanList(nc)
                        % make plot
                        inputData.filtered = 1;
                        inputData.plotAllChannels = 0;
                        plotArtifacts_singleParameter(inputData,artifactData,j,stimChanList(nc),nw);
                    end
                end
            end
        end
    else % make plots for each artifactData, assume cathodal and anodal interleaved
        % make a 10x10 figure with the response for each response channel
        for i = 1:numel(artifactData)
            outputFigures(end+1)=figure;
            set(outputFigures(end),'Name',['StimOnCH-',num2str(artifactData(i).stimChannel),'_AllChannels_FileNum-',num2str(i)])
            numPlotPixels=1200;
            set(outputFigures(end),'Position',[100 100 numPlotPixels numPlotPixels]);
            paperSize=0.2+numPlotPixels/get(outputFigures(end),'ScreenPixelsPerInch');
            set(outputFigures(end),'PaperSize',[paperSize,paperSize]);
            for j = 1:numel(artifactData(i).electrodeNames)
                %find the index in our full list of electrodes, that
                %corresponds to the channel j in artifactData(i)
                posIdx=find(strcmp(eList,artifactData(i).electrodeNames{j}));
                eRow=posList(posIdx,1);
                eCol=posList(posIdx,2);
                h=subplot(10,10,10*(eRow-1)+eCol);
                hold on
                if chList(posIdx)==artifactData(i).stimChannel
                    %put a purple box around the stim channel:
                    stimBoxX=[2,.95*size(artifactData(i).artifact,3),.95*size(artifactData(i).artifact,3),2,2];
                    stimBoxY=[-inputData.plotRange*.95*1000,-inputData.plotRange*.95*1000,inputData.plotRange*.95*1000,inputData.plotRange*.95*1000,-inputData.plotRange*.95*1000];
                    plot(stimBoxX,stimBoxY,'mp-','linewidth',3)
                end
                if find(chList(posIdx)==inputData.badChList,1,'first')
                   %put a red X through known bad channels:
                   badMarkX=[2,.95*size(artifactData(i).artifact,3)];
                   badMarkY=[-inputData.plotRange*.95*1000,inputData.plotRange*.95*1000];
                   plot(badMarkX,badMarkY,'rp-','lineWidth',3)
                   badMarkX=[2,.95*size(artifactData(i).artifact,3)];
                   badMarkY=[inputData.plotRange*.95*1000,-inputData.plotRange*.95*1000];
                   plot(badMarkX,badMarkY,'rp-','lineWidth',3)
                end
                inputData.filtered = 0;
                inputData.plotAllChannels = 1;
                plotArtifacts_cathodalAnodal(inputData,artifactData,i,j);
            end
        end
        
        % make a single figure focusing on the stimulated channel
        for i = 1:numel(artifactData)
            outputFigures(end+1)=figure;
            set(outputFigures(end),'Name',['StimOnCH-',num2str(artifactData(i).stimChannel),'_Raw_FileNum-',num2str(i)])
            for j = 1:numel(artifactData(i).electrodeNames)
                posIdx=find(strcmp(eList,artifactData(i).electrodeNames{j}));
                if(chList(posIdx)==artifactData(i).stimChannel)
                    inputData.filtered = 0;
                    inputData.plotAllChannels = 0;
                    plotArtifacts_cathodalAnodal(inputData,artifactData,i,j);
                end
            end
        end
        
        % make a single figure focusing on the stimulated channel --
        % filtered
        for i = 1:numel(artifactData)
            outputFigures(end+1)=figure;
            set(outputFigures(end),'Name',['StimOnCH-',num2str(artifactData(i).stimChannel),'_Filtered_FileNum-',num2str(i)])
            for j = 1:numel(artifactData(i).electrodeNames)
                posIdx=find(strcmp(eList,artifactData(i).electrodeNames{j}));
                if chList(posIdx)==artifactData(i).stimChannel
                    % plot data
                    inputData.filtered = 1;
                    inputData.plotAllChannels = 0;
                    plotArtifacts_cathodalAnodal(inputData,artifactData,i,j);
                end
            end
        end
    end

    if(inputData.saveFigures == 0)
        outputFigures = [];
    end
end

function [] = plotArtifacts_cathodalAnodal(inputData,artifactData,artifactDataIdx,chanIdx)
    %identify which artifacts are from cathodal first pulses:
    cathodalMask=false(numel(artifactData(artifactDataIdx).stimOn),1);
    cathodalMask([1:10:min(100,numel(artifactData(artifactDataIdx).stimOn))])=true;
    anodalMask=false(numel(artifactData(artifactDataIdx).stimOn),1);
    anodalMask([2:10:min(100,numel(artifactData(artifactDataIdx).stimOn))])=true;
    %make colormasks for the  cathodal and anodal pulses

    %plot the pulses: cathodal red, anodal blue
    oneFill=ones(sum(cathodalMask),1);
    halfFall=.75-.75/sum(cathodalMask):-.75/sum(cathodalMask):0;
    stimCathColor=[oneFill,halfFall',halfFall'];
    set(gca,'colororder',stimCathColor);
    hold all
    xData = ((1:size(artifactData(artifactDataIdx).artifact,3))-1)/30 - inputData.presample/30;
    if(inputData.filtered)
        plot(xData,squeeze(artifactData(artifactDataIdx).artifactFiltered(chanIdx,cathodalMask,:))')
    else
        plot(xData,squeeze(artifactData(artifactDataIdx).artifact(chanIdx,cathodalMask,:))'-squeeze(repmat(artifactData(artifactDataIdx).artifact(chanIdx,cathodalMask,1),[1,1,size(artifactData(artifactDataIdx).artifact,3)]))')
    end
    %set colororder to blue gradient and plot anodal pulses:
    halfFall=.5-.5/sum(~cathodalMask):-.5/sum(~cathodalMask):0;
    oneFill=ones(sum(~cathodalMask),1);
    stimAnColor=[halfFall',halfFall',oneFill];
    set(gca,'colororder',stimAnColor);
    hold all
    
    if(inputData.filtered)
        plot(xData,squeeze(artifactData(artifactDataIdx).artifactFiltered(chanIdx,anodalMask,:))')
    else
        plot(xData,squeeze(artifactData(artifactDataIdx).artifact(chanIdx,anodalMask,:))'-squeeze(repmat(artifactData(artifactDataIdx).artifact(chanIdx,anodalMask,1),[1,1,size(artifactData(artifactDataIdx).artifact,3)]))')
    end
    
    axis tight%keeps from padding time, we will set y axis below:
    if(inputData.filtered)
        ylim([-inputData.plotRangeFiltered,inputData.plotRangeFiltered]*1000);
    else
        ylim([-inputData.plotRange,inputData.plotRange]*1000);
    end
    
    if(inputData.plotAllChannels)
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
    else
        formatForLee(gcf)
    end
end     
             

function [] = plotArtifacts_singleParameter(inputData,artifactData,chanIdx,chanNum,waveNum)
    %identify which artifacts are from the channel and wave num:
    stimMask=false(numel(artifactData(1).stimOn),1);
    stimMask(artifactData(1).chanSent == chanNum & artifactData(1).waveNum == waveNum)=true;
    if(sum(stimMask) > 10)
        stimMask(find(stimMask,sum(stimMask)-10,'last')) = 0;
    end
    %make colormasks for the pulses
    %plot the pulses: black
    oneFill=ones(sum(stimMask),1);
    halfFall=.75-.75/sum(stimMask):-.75/sum(stimMask):0;
    stimColor=[halfFall',halfFall',halfFall'];
    set(gca,'colororder',stimColor);
    hold all
    
    xData = ((1:size(artifactData(1).artifact,3))-1)/30 - inputData.presample/30;
    if(inputData.filtered)
        plot(xData,squeeze(artifactData(1).artifactFiltered(chanIdx,stimMask,:))')
    else
        plot(xData,squeeze(artifactData(1).artifact(chanIdx,stimMask,:))'-squeeze(repmat(artifactData(1).artifact(chanIdx,stimMask,1),[1,1,size(artifactData(1).artifact,3)]))')
    end
    
    axis tight%keeps from padding time, we will set y axis below:
    if(inputData.filtered)
        ylim([-inputData.plotRangeFiltered,inputData.plotRangeFiltered]*1000);
    else
        ylim([-inputData.plotRange,inputData.plotRange]*1000);
    end
    if(inputData.plotAllChannels)
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
    else
        formatForLee(gcf)
    end

end
