function [ outputFigures,outputData ] = processStimArtifact(folderpath, inputData )
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
    
%     %% plot the artifacts for each stim channel:
% %     stimChansFig=figure;
%     for i=1:numel(artifactData)
%     %    disp(['plotting stim:',num2str(i)])
%         %loop across the data for each response channel and put it into a
%         %subplot:
%         outputFigures(end+1)=figure;
%         set(outputFigures(end),'Name',['StimOnCH-',num2str(artifactData(i).stimChannel)])
%         numPlotPixels=1200;
%         set(outputFigures(end),'Position',[100 100 numPlotPixels numPlotPixels]);
%         paperSize=0.2+numPlotPixels/get(outputFigures(end),'ScreenPixelsPerInch');
%         set(outputFigures(end),'PaperSize',[paperSize,paperSize]);
%         for j=1:numel(artifactData(i).electrodeNames)
%             %disp(['   plotting response for electrode:',num2str(j)])
%             %find the index in our full list of electrodes, that
%             %corresponds to the channel j in artifactData(i)
%             posIdx=find(strcmp(eList,artifactData(i).electrodeNames{j}));
%             eRow=posList(posIdx,1);
%             eCol=posList(posIdx,2);
%             h=subplot(10,10,10*(eRow-1)+eCol);
%             hold on
%             %add background fill for stimulus:
%             %assumes 200us/phase+53us interpulse and rounds up to the next
%             %1s/30000 tick
%             stimLength=ceil((inputData.pWidth1+inputData.pWidth2+inputData.interpulse)*30000);
%             stimFillBox=[inputData.presample, -8200;...
%                             inputData.presample, 8200;...
%                             inputData.presample+stimLength, 8200;...
%                             inputData.presample+stimLength, -8200;...
%                             inputData.presample, -8200];
%             stimFillColor=[1 .75 .75];
%             fill(stimFillBox(:,1),stimFillBox(:,2),stimFillColor);
%             %find the sync-off time and put a fill box out to the sync off
%             %time:
%             stimOffOffset=ceil(nanmean(artifactData(i).stimOff-artifactData(i).stimOn));
%             %add background fill for post-stim sync on period
%             syncFillBox=[inputData.presample+stimLength,-8200;...
%                             inputData.presample+stimOffOffset,-8200;...
%                             inputData.presample+stimOffOffset,8200;...
%                             inputData.presample+stimLength,8200;...
%                             inputData.presample+stimLength,-8200];
%             syncFillColor=[1 1 .5];
%             fill(syncFillBox(:,1),syncFillBox(:,2),syncFillColor)
%             %add line @ 1ms:
%             plot([1,1]*(inputData.presample+30),[8200,-8200],'g')
%             %add line @ 1.5ms:
%             plot([1,1]*(inputData.presample+45),[8200,-8200],'b')
%             %identify which artifacts are from cathodal first pulses:
%             cathodalMask=false(numel(artifactData(i).stimOn),1);
%             cathodalMask([1:10:numel(artifactData(i).stimOn)])=true;
%             anodalMask=false(numel(artifactData(i).stimOn),1);
%             anodalMask([2:10:numel(artifactData(i).stimOn)])=true;
%             %make corlormasks for the  cathodal and anodal pulses
%             if chList(posIdx)==artifactData(i).stimChannel
%                 %put a purple box around the stim channel:
%                 stimBoxX=[2,.95*size(artifactData(i).artifact,3),.95*size(artifactData(i).artifact,3),2,2];
%                 stimBoxY=[-inputData.plotRange*.95*1000,-inputData.plotRange*.95*1000,inputData.plotRange*.95*1000,inputData.plotRange*.95*1000,-inputData.plotRange*.95*1000];
%                 plot(stimBoxX,stimBoxY,'mp-','linewidth',3)
%             end
%             if find(chList(posIdx)==inputData.badChList,1,'first')
%                %put a red X through known bad channels:
%                badMarkX=[2,.95*size(artifactData(i).artifact,3)];
%                badMarkY=[-inputData.plotRange*.95*1000,inputData.plotRange*.95*1000];
%                plot(badMarkX,badMarkY,'rp-','lineWidth',3)
%                badMarkX=[2,.95*size(artifactData(i).artifact,3)];
%                badMarkY=[inputData.plotRange*.95*1000,-inputData.plotRange*.95*1000];
%                plot(badMarkX,badMarkY,'rp-','lineWidth',3)
%             end
%             %plot the pulses: cathodal red, anodal blue
%             oneFill=ones(sum(cathodalMask),1);
%             halfFall=.75-.75/sum(cathodalMask):-.75/sum(cathodalMask):0;
%             stimCathColor=[oneFill,halfFall',halfFall'];
%             set(gca,'colororder',stimCathColor);
%             hold all
%             plot(squeeze(artifactData(i).artifact(j,cathodalMask,:))'-squeeze(repmat(artifactData(i).artifact(j,cathodalMask,1),[1,1,size(artifactData(i).artifact,3)]))')
%             %set colororder to blue gradient and plot anodal pulses:
%             halfFall=.5-.5/sum(~cathodalMask):-.5/sum(~cathodalMask):0;
%             oneFill=ones(sum(~cathodalMask),1);
%             stimAnColor=[halfFall',halfFall',oneFill];
%             set(gca,'colororder',stimAnColor);
%             hold all
%             plot(squeeze(artifactData(i).artifact(j,anodalMask,:))'-squeeze(repmat(artifactData(i).artifact(j,anodalMask,1),[1,1,size(artifactData(i).artifact,3)]))')
%             axis tight%keeps from padding time, we will set y axis below:
%             ylim([-inputData.plotRange,inputData.plotRange]*1000);
%             set(gca,'XTickLabel',[])
%             set(gca,'YTickLabel',[])
%             
%         end
%         set(gcf,'NextPlot','add');
%         axes;
%         h = title(['Artifact in response to stim on Channel: ',num2str(artifactData(i).stimChannel)]);
%         set(gca,'Visible','off');
%         set(h,'Visible','on'); 
%     end
%             
%     %% now run through and plot the monitor data if it exists
%     if exist('mIdx','var')
%         for i=1:numel(artifactData)
%             outputFigures(end+1)=figure;
%             set(outputFigures(end),'Name',['MonitorCH-',num2str(artifactData(i).stimChannel)])
%             numPlotPixels=1200;
%             set(outputFigures(end),'Position',[100 100 numPlotPixels numPlotPixels]);
%             paperSize=0.2+numPlotPixels/get(outputFigures(end),'ScreenPixelsPerInch');
%             set(outputFigures(end),'PaperSize',[paperSize,paperSize]);
%             for j=1:numel(artifactData(i).electrodeNames)
%                 %disp(['   plotting response for electrode:',num2str(j)])
%                 %find the index in our full list of electrodes, that
%                 %corresponds to the channel j in artifactData(i)
%                 posIdx=find(strcmp(eList,artifactData(i).electrodeNames{j}));
%                 eRow=posList(posIdx,1);
%                 eCol=posList(posIdx,2);
%                 h=subplot(10,10,10*(eRow-1)+eCol);
%                 hold on
%                 %add background fill for stimulus:
%                 %assumes 200us/phase+53us interpulse and rounds up to the next
%                 %1s/30000 tick
%                 stimLength=ceil((.0002*2+.000053)*30000);
%                 stimFillBox=[inputData.presample, -8200;...
%                                 inputData.presample, 8200;...
%                                 inputData.presample+stimLength, 8200;...
%                                 inputData.presample+stimLength, -8200;...
%                                 inputData.presample, -8200];
%                 stimFillColor=[1 .75 .75];
%                 fill(stimFillBox(:,1),stimFillBox(:,2),stimFillColor);
%                 %find the sync-off time and put a fill box out to the sync off
%                 %time:
%                 stimOffOffset=ceil(nanmean(artifactData(i).stimOff-artifactData(i).stimOn));
%                 %add background fill for post-stim sync on period
%                 syncFillBox=[inputData.presample+stimLength,-8200;...
%                                 inputData.presample+stimOffOffset,-8200;...
%                                 inputData.presample+stimOffOffset,8200;...
%                                 inputData.presample+stimLength,8200;...
%                                 inputData.presample+stimLength,-8200];
%                 syncFillColor=[1 1 .5];
%                 fill(syncFillBox(:,1),syncFillBox(:,2),syncFillColor)
%                 %add line @ 1ms:
%                 plot([1,1]*(inputData.presample+30),[8200,-8200],'g')
%                 %add line @ 1.5ms:
%                 plot([1,1]*(inputData.presample+45),[8200,-8200],'b')
%                 %identify which monitor channel data is from cathodal first pulses:
%                 range1=mean(max(artifactData(i).monitor1(1:2:end),[],2)-min(artifactData(i).monitor1(1:2:end),[],2));
%                 range2=mean(max(artifactData(i).monitor1(2:2:end),[],2)-min(artifactData(i).monitor1(2:2:end),[],2));
%                 if range1>range2
%                     %monitor1 is the cathodal pulse:
%                     cathodalPulses=squeeze(artifactData(i).monitor1(1:10:end,:));
%                     anodalPulses=squeeze(artifactData(i).monitor2(2:10:end,:));
%                 else
%                    %monitor1 is the anodal pulse:
%                     cathodalPulses=squeeze(artifactData(i).monitor2(1:10:end,:));
%                     anodalPulses=squeeze(artifactData(i).monitor1(2:10:end,:));
%                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                 cathodalMask=false(numel(artifactData(i).stimOn),1);
%                 cathodalMask([1:10:numel(artifactData(i).stimOn)])=true;
%                 anodalMask=false(numel(artifactData(i).stimOn),1);
%                 anodalMask([2:10:numel(artifactData(i).stimOn)])=true;
%                 %make corlormasks for the  cathodal and anodal pulses
%                 if chList(posIdx)==artifactData(i).stimChannel
%                     %put a purple box around the stim channel:
%                     stimBoxX=[2,                    .95*size(artifactData(i).artifact,3),    .95*size(artifactData(i).artifact,3),    2,                      2];
%                     stimBoxY=[-inputData.plotRange*.95*1000, -inputData.plotRange*.95*1000,       inputData.plotRange*.95*1000,        inputData.plotRange*.95*1000,    -inputData.plotRange*.95*1000];
%                     plot(stimBoxX,stimBoxY,'mp-','linewidth',3)
%                 end
%                 %plot the pulses: cathodal red, anodal blue
%                 oneFill=ones(sum(cathodalMask),1);
%                 halfFall=.75-.75/sum(cathodalMask):-.75/sum(cathodalMask):0;
%                 stimCathColor=[oneFill,halfFall',halfFall'];
%                 set(gca,'colororder',stimCathColor);
%                 hold all
%                 plot(squeeze(artifactData(i).artifact(j,cathodalMask,:))'-squeeze(repmat(artifactData(i).artifact(j,cathodalMask,1),[1,1,size(artifactData(i).artifact,3)]))')
%                 %set colororder to blue gradient and plot anodal pulses:
%                 halfFall=.5-.5/sum(~cathodalMask):-.5/sum(~cathodalMask):0;
%                 oneFill=ones(sum(~cathodalMask),1);
%                 stimAnColor=[halfFall',halfFall',oneFill];
%                 set(gca,'colororder',stimAnColor);
%                 hold all
%                 plot(plot(squeeze(artifactData(i).artifact(j,anodalMask,:))'-squeeze(repmat(artifactData(i).artifact(j,anodalMask,1),[1,1,size(artifactData(i).artifact,3)]))'))
%                 axis tight%keeps from padding time, we will set y axis below:
%                 ylim([-inputData.plotRange,inputData.plotRange]*1000);
%                 set(gca,'XTickLabel',[])
%                 set(gca,'YTickLabel',[])
% 
%             end
%             set(gcf,'NextPlot','add');
%             axes;
%             h = title(['Artifact in response to stim on Channel: ',num2str(artifactData(i).stimChannel)]);
%             set(gca,'Visible','off');
%             set(h,'Visible','on'); 
%         end
%     end
%     %% now export figures directly to png so we don't hit the print command quite so often
%     for i=1:numel(outputFigures)
%         %save a png of the figure
%         fname=get(outputFigures(i),'Name');
%         if isempty(fname)
%             fname=strcat('Figure_',num2str(double(H)));
%         end
%         fname(fname==' ')='_';%replace spaces in name for saving
%         print('-dpng',outputFigures(i),strcat(folderpath,['Raw_Figures' filesep 'PNG' filesep],fname,'.png'))
%     end
    outputFigures=[];
end

