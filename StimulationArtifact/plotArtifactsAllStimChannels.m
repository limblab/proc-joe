<<<<<<< HEAD
function [  ] = plotArtifactsAllStimChannels( outputData, inputData, folderpath, varargin )

figName = '';
noPlots = 0;
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'Name'
            figName = varargin{i+1};
        case 'noPlots'
            noPlots = varargin{i+1};
    end
end

outputFigures = [];
artifactData = outputData.artifactData;
eList = outputData.eList;
posList = outputData.posList;
chList = outputData.chList;

%% plot the artifacts for each stim channel:
%     stimChansFig=figure;
    for i=1:numel(artifactData)
    %    disp(['plotting stim:',num2str(i)])
        %loop across the data for each response channel and put it into a
        %subplot:
        if(noPlots)
            outputFigures(end+1)=figure('Visible','off');
        else
            outputFigures(end+1)=figure('Visible','on');
        end
        set(outputFigures(end),'Name',['StimOnCH-',num2str(artifactData(i).stimChannel),figName])
        numPlotPixels=1200;
        set(outputFigures(end),'Position',[100 100 numPlotPixels numPlotPixels]);
        paperSize=0.2+numPlotPixels/get(outputFigures(end),'ScreenPixelsPerInch');
        set(outputFigures(end),'PaperSize',[paperSize,paperSize]);
        for j=1:numel(artifactData(i).electrodeNames)
            %disp(['   plotting response for electrode:',num2str(j)])
            %find the index in our full list of electrodes, that
            %corresponds to the channel j in artifactData(i)
            posIdx=find(strcmp(eList,artifactData(i).electrodeNames{j}));
            eRow=posList(posIdx,1);
            eCol=posList(posIdx,2);
            h=subplot(10,10,10*(eRow-1)+eCol);
            hold on
            %add background fill for stimulus:
            %assumes 200us/phase+53us interpulse and rounds up to the next
            %1s/30000 tick
            stimLength=ceil((inputData.pWidth1+inputData.pWidth2+inputData.interpulse)*30000);
            stimFillBox=[inputData.presample, -8200;...
                            inputData.presample, 8200;...
                            inputData.presample+stimLength, 8200;...
                            inputData.presample+stimLength, -8200;...
                            inputData.presample, -8200];
            stimFillColor=[1 .75 .75];
            fill(stimFillBox(:,1),stimFillBox(:,2),stimFillColor);
            %find the sync-off time and put a fill box out to the sync off
            %time:
            stimOffOffset=ceil(nanmean(artifactData(i).stimOff-artifactData(i).stimOn));
            %add background fill for post-stim sync on period
            syncFillBox=[inputData.presample+stimLength,-8200;...
                            inputData.presample+stimOffOffset,-8200;...
                            inputData.presample+stimOffOffset,8200;...
                            inputData.presample+stimLength,8200;...
                            inputData.presample+stimLength,-8200];
            syncFillColor=[1 1 .5];
            fill(syncFillBox(:,1),syncFillBox(:,2),syncFillColor)
            %add line @ 1ms:
            plot([1,1]*(inputData.presample+30),[8200,-8200],'g')
            %add line @ 1.5ms:
            plot([1,1]*(inputData.presample+45),[8200,-8200],'b')
            %identify which artifacts are from cathodal first pulses:
            cathodalMask=false(numel(artifactData(i).stimOn),1);
            cathodalMask([1:10:numel(artifactData(i).stimOn)])=true;
            anodalMask=false(numel(artifactData(i).stimOn),1);
            anodalMask([2:10:numel(artifactData(i).stimOn)])=true;
            %make corlormasks for the  cathodal and anodal pulses
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
            %plot the pulses: cathodal red, anodal blue
            oneFill=ones(sum(cathodalMask),1);
            halfFall=.75-.75/sum(cathodalMask):-.75/sum(cathodalMask):0;
            stimCathColor=[oneFill,halfFall',halfFall'];
            set(gca,'colororder',stimCathColor);
            hold all
            plot(squeeze(artifactData(i).artifact(j,cathodalMask,:))'-squeeze(repmat(artifactData(i).artifact(j,cathodalMask,1),[1,1,size(artifactData(i).artifact,3)]))')
            %set colororder to blue gradient and plot anodal pulses:
            halfFall=.5-.5/sum(~cathodalMask):-.5/sum(~cathodalMask):0;
            oneFill=ones(sum(~cathodalMask),1);
            stimAnColor=[halfFall',halfFall',oneFill];
            set(gca,'colororder',stimAnColor);
            hold all
            plot(squeeze(artifactData(i).artifact(j,anodalMask,:))'-squeeze(repmat(artifactData(i).artifact(j,anodalMask,1),[1,1,size(artifactData(i).artifact,3)]))')
            axis tight%keeps from padding time, we will set y axis below:
            ylim([-inputData.plotRange,inputData.plotRange]*1000);
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
            
        end
        set(gcf,'NextPlot','add');
        axes;
        figTitle = get(outputFigures(end),'Name');
        figTitle = strrep(figTitle,'_','\_');% add \ before _
        
        h = title(figTitle);
        set(gca,'Visible','off');
        set(h,'Visible','on'); 
    end
    %% now export figures directly to png so we don't hit the print command quite so often
    for i=1:numel(outputFigures)
        %save a png of the figure
        fname=get(outputFigures(i),'Name');
        if isempty(fname)
            fname=strcat('Figure_',num2str(double(H)));
        end
        fname(fname==' ')='_';%replace spaces in name for saving
        print('-dpng',outputFigures(i),strcat(cd,filesep,'Raw_Figures',filesep,fname,'.png'))
    end
    outputFigures=[];

end

=======
function [  ] = plotArtifactsAllStimChannels( outputData, inputData, folderpath, varargin )

figName = '';
noPlots = 0;
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'Name'
            figName = varargin{i+1};
        case 'noPlots'
            noPlots = varargin{i+1};
    end
end

outputFigures = [];
artifactData = outputData.artifactData;
eList = outputData.eList;
posList = outputData.posList;
chList = outputData.chList;

%% plot the artifacts for each stim channel:
%     stimChansFig=figure;
    for i=1:numel(artifactData)
    %    disp(['plotting stim:',num2str(i)])
        %loop across the data for each response channel and put it into a
        %subplot:
        if(noPlots)
            outputFigures(end+1)=figure('Visible','off');
        else
            outputFigures(end+1)=figure('Visible','on');
        end
        set(outputFigures(end),'Name',['StimOnCH-',num2str(artifactData(i).stimChannel),figName])
        numPlotPixels=1200;
        set(outputFigures(end),'Position',[100 100 numPlotPixels numPlotPixels]);
        paperSize=0.2+numPlotPixels/get(outputFigures(end),'ScreenPixelsPerInch');
        set(outputFigures(end),'PaperSize',[paperSize,paperSize]);
        for j=1:numel(artifactData(i).electrodeNames)
            %disp(['   plotting response for electrode:',num2str(j)])
            %find the index in our full list of electrodes, that
            %corresponds to the channel j in artifactData(i)
            posIdx=find(strcmp(eList,artifactData(i).electrodeNames{j}));
            eRow=posList(posIdx,1);
            eCol=posList(posIdx,2);
            h=subplot(10,10,10*(eRow-1)+eCol);
            hold on
            %add background fill for stimulus:
            %assumes 200us/phase+53us interpulse and rounds up to the next
            %1s/30000 tick
            stimLength=ceil((inputData.pWidth1+inputData.pWidth2+inputData.interpulse)*30000);
            stimFillBox=[inputData.presample, -8200;...
                            inputData.presample, 8200;...
                            inputData.presample+stimLength, 8200;...
                            inputData.presample+stimLength, -8200;...
                            inputData.presample, -8200];
            stimFillColor=[1 .75 .75];
            fill(stimFillBox(:,1),stimFillBox(:,2),stimFillColor);
            %find the sync-off time and put a fill box out to the sync off
            %time:
            stimOffOffset=ceil(nanmean(artifactData(i).stimOff-artifactData(i).stimOn));
            %add background fill for post-stim sync on period
            syncFillBox=[inputData.presample+stimLength,-8200;...
                            inputData.presample+stimOffOffset,-8200;...
                            inputData.presample+stimOffOffset,8200;...
                            inputData.presample+stimLength,8200;...
                            inputData.presample+stimLength,-8200];
            syncFillColor=[1 1 .5];
            fill(syncFillBox(:,1),syncFillBox(:,2),syncFillColor)
            %add line @ 1ms:
            plot([1,1]*(inputData.presample+30),[8200,-8200],'g')
            %add line @ 1.5ms:
            plot([1,1]*(inputData.presample+45),[8200,-8200],'b')
            %identify which artifacts are from cathodal first pulses:
            cathodalMask=false(numel(artifactData(i).stimOn),1);
            cathodalMask([1:10:numel(artifactData(i).stimOn)])=true;
            anodalMask=false(numel(artifactData(i).stimOn),1);
            anodalMask([2:10:numel(artifactData(i).stimOn)])=true;
            %make corlormasks for the  cathodal and anodal pulses
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
            %plot the pulses: cathodal red, anodal blue
            oneFill=ones(sum(cathodalMask),1);
            halfFall=.75-.75/sum(cathodalMask):-.75/sum(cathodalMask):0;
            stimCathColor=[oneFill,halfFall',halfFall'];
            set(gca,'colororder',stimCathColor);
            hold all
            plot(squeeze(artifactData(i).artifact(j,cathodalMask,:))'-squeeze(repmat(artifactData(i).artifact(j,cathodalMask,1),[1,1,size(artifactData(i).artifact,3)]))')
            %set colororder to blue gradient and plot anodal pulses:
            halfFall=.5-.5/sum(~cathodalMask):-.5/sum(~cathodalMask):0;
            oneFill=ones(sum(~cathodalMask),1);
            stimAnColor=[halfFall',halfFall',oneFill];
            set(gca,'colororder',stimAnColor);
            hold all
            plot(squeeze(artifactData(i).artifact(j,anodalMask,:))'-squeeze(repmat(artifactData(i).artifact(j,anodalMask,1),[1,1,size(artifactData(i).artifact,3)]))')
            axis tight%keeps from padding time, we will set y axis below:
            ylim([-inputData.plotRange,inputData.plotRange]*1000);
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
            
        end
        set(gcf,'NextPlot','add');
        axes;
        figTitle = get(outputFigures(end),'Name');
        figTitle = strrep(figTitle,'_','\_');% add \ before _
        
        h = title(figTitle);
        set(gca,'Visible','off');
        set(h,'Visible','on'); 
    end
    %% now export figures directly to png so we don't hit the print command quite so often
    for i=1:numel(outputFigures)
        %save a png of the figure
        fname=get(outputFigures(i),'Name');
        if isempty(fname)
            fname=strcat('Figure_',num2str(double(H)));
        end
        fname(fname==' ')='_';%replace spaces in name for saving
        print('-dpng',outputFigures(i),strcat(folderpath,filesep,'Raw_Figures',filesep,fname,'.png'))
    end
    outputFigures=[];

end

>>>>>>> 0e441d3b4e9a2e8831546d9f89425fae8b961b8d
