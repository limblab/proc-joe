% plot artifacts on a subset of the array

figure('Position',[357.8000 56.2000 838.4000 726.4000]);
mapData=loadMapFile(mapFileName);
%establish trackers for information that is common to all files:
posList=[];
posList = [mapData.row,mapData.col];
plottedHere = zeros(10,10);% plot neurons
stimChan = 73;

for chan = 1:96
    posIdx=find(mapData.chan==chan);
    eRow=posList(posIdx,1);
    eRow = 11 - eRow;
    eCol=posList(posIdx,2);
%     h=subplot(10,10,10*(eRow-1)+eCol);
%     plottedHere(eRow,eCol) = 1;

    if(eRow >= 4 && eRow <= 8 && eCol <= 3) 
        h=subplot(5,3,3*(eRow-4)+eCol);
        plottedHere(eRow,eCol) = 1;
        %% PLOT ARTIFACTS CODE
        xData = ((0:1:size(dataStruct2.artifactData(1).artifact,3)-1)-inputData.presample)/30;
        plot(xData,squeeze(dataStruct2.artifactData(1).artifact(chan,1:10:end,:))','r','linewidth',1.5)
        hold on
        plot(xData,squeeze(dataStruct2.artifactData(1).artifact(chan,2:10:end,:))','b','linewidth',1.5)

        % remove axes
        % plot one set of axes
        xlim([-0.3,5])
        ylim([-1000,1000])
        if(eRow == 8 && eCol == 1)
            set(gca,'visible','on')
            set(gca,'fontsize',16)
            formatForLee(gcf)
            set(gca,'YTick',[-1000,1000])
            set(gca,'YMinorTick','off')
            set(gca,'YTickLabel',[])
            xlabel('Time after stimulation onset (ms)')
            ylabel('Voltage (\muV)')
        else
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
            set(gca,'XLabel',[]);
            set(gca,'YLabel',[]);
            set(gca,'XTick',[]);
            set(gca,'XMinorTick','off')
            set(gca,'YTick',[]);
            set(gca,'YMinorTick','off')
            set(gca,'Visible','off')
        end
        
        
        
    end
end

% mark stim chan
posIdx=find(mapData.chan==stimChan);
eRow=posList(posIdx,1);
eRow = 11 - eRow;
eCol=posList(posIdx,2);
% h=subplot(10,10,10*(eRow-1)+eCol);
h=subplot(5,3,3*(eRow-4)+eCol);
hold on
ax = gca;
XLIM = ax.XLim;
YLIM = ax.YLim;
if(plottedHere(eRow,eCol)==1)
    stimBoxX=[XLIM(1),XLIM(2),XLIM(2),XLIM(1),XLIM(1)];
    stimBoxY=[YLIM(1),YLIM(1),YLIM(2),YLIM(2),YLIM(1)];
    plot(stimBoxX,stimBoxY,'m-','linewidth',3)
else
%                 stimCrossX1 = [XLIM(1),XLIM(2)];
%                 stimCrossY1 = [YLIM(1),YLIM(2)];
%                 stimCrossX2 = [XLIM(1),XLIM(2)];
%                 stimCrossY2 = [YLIM(2),YLIM(1)];
    stimDotX = (XLIM(2) + XLIM(1))/2;
    stimDotY = (YLIM(2) + YLIM(1))/2;
    plot(stimDotX,stimDotY,'m','markersize',75)
    xlim(XLIM)
    ylim(YLIM)
    set(ax,'visible','off')
    set(ax,'color','none')
end
