function [ output_args ] = plotPSTHStimWholeArray( cds,mapFilename, varargin )
% plot PSTH for the whole array

plotLine = 1;
preTime = 10/1000;
postTime = 30/1000;
binSize = 0.0002;
waveformTypes = 1:numel(cds.waveforms.parameters);
chans = 1:numel(unique(cds.waveforms.chanSent));
lineColor = 'k';
lineWidth = 1;
plotAllOnOneFigure = 0;
plotStimOn = 1;
plotStimChan = 1;
plotProbabilityText = 0;
for i = 1:2:numel(varargin)
    switch varargin{i}
        case 'plotLine'
            plotLine = varargin{i+1};
        case 'lineColor'
            lineColor = varargin{i+1};
        case 'plotAllOnOneFigure'
            plotAllOnOneFigure = varargin{i+1};
        case 'preTime'
            preTime = varargin{i+1};
        case 'postTime'
            postTime = varargin{i+1};
        case 'binSize'
            binSize = varargin{i+1};
        case 'waveformTypes'
            waveformTypes = varargin{i+1};
        case 'chans'
            chans = varargin{i+1};
        case 'plotStimOn'
            plotStimOn = varargin{i+1};
        case 'plotStimChan'
            plotStimChan = varargin{i+1};
        case 'plotProbabilityText'
            plotProbabilityText = varargin{i+1};
    end
end

mapData=loadMapFile(mapFilename);

%establish trackers for information that is common to all files:
posList=[];
posList = [mapData.row,mapData.col];
plottedHere = zeros(10,10);
for c = 1:numel(chans)
    for wave = 1:numel(waveformTypes)
        figure();
        % plot neurons
        for nn = 1:size(cds.units,2)
            if(cds.units(nn).ID~=0 && cds.units(nn).ID~=255)
                posIdx=find(mapData.chan==cds.units(nn).chan);
                eRow=posList(posIdx,1);
                eRow = 11 - eRow;
                eCol=posList(posIdx,2);
                h=subplot(10,10,10*(eRow-1)+eCol);
                plottedHere(eRow,eCol) = 1;
                
                plotPSTHStim(cds,nn,'binSize',binSize,'makeFigure',0,'makeSubplots',0,'plotTitle',0,'waveformTypes',waveformTypes(wave),...
                        'chans',chans(c),'preTime',preTime,'postTime',postTime,'lineColor',lineColor,'plotLine',plotLine,'plotStimOn',plotStimOn,...
                        'lineWidth',lineWidth);

                % remove axes
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
                set(gca,'XLabel',[]);
                set(gca,'YLabel',[]);
                set(gca,'XTick',[]);
                set(gca,'XMinorTick','off')
                set(gca,'YTick',[]);
                set(gca,'YMinorTick','off')
                if(plotProbabilityText)
                    prob = getProbabilityOfResponse(cds,nn,'peakPeriod','automatic','preTime',20/1000,'postTime',60/1000,...
                        'chans',chans(c),'waveformTypes',waveformTypes(wave));
                    if(prob > 0)
                        
                        ax = gca;
                        XLIM = ax.XLim;
                        YLIM = ax.YLim;
                        text(XLIM(2)*0.75,YLIM(2)*1.25,strcat(num2str(round(prob*100,1)),'%'),'color','b')
                    end
                end
            end
        end
        
        % plot X through channels without neurons
        for i = 1:size(plottedHere,1)
            for j = 1:size(plottedHere,2)
                if(plottedHere(i,j) == 0)
                    h=subplot(10,10,10*(i-1)+j);
                    hold on
                    ax = gca;
                    XLIM = ax.XLim;
                    YLIM = ax.YLim;
%                     stimCrossX1 = [XLIM(1),XLIM(2)];
%                     stimCrossY1 = [YLIM(1),YLIM(2)];
%                     stimCrossX2 = [XLIM(1),XLIM(2)];
%                     stimCrossY2 = [YLIM(2),YLIM(1)];
%                     plot(stimCrossX1,stimCrossY1,'kp-','linewidth',3)
%                     plot(stimCrossX2,stimCrossY2,'kp-','linewidth',3)
                    stimDotX = (XLIM(2) + XLIM(1))/2;
                    stimDotY = (YLIM(2) + YLIM(1))/2;
                    plot(stimDotX,stimDotY,'k.','markersize',75)
                    xlim(XLIM)
                    ylim(YLIM)
                    set(ax,'visible','off')
                    set(ax,'color','none')
                end
            end
        end
        
        if(plotStimChan)
            chanSent = unique(cds.waveforms.chanSent);
            stimChan = chanSent(chans(c));
            posIdx=find(mapData.chan==stimChan);
            eRow=posList(posIdx,1);
            eRow = 11 - eRow;
            eCol=posList(posIdx,2);
            h=subplot(10,10,10*(eRow-1)+eCol);
            hold on
            ax = gca;
            XLIM = ax.XLim;
            YLIM = ax.YLim;
            if(plottedHere(eRow,eCol)==1)
                stimBoxX=[XLIM(1),XLIM(2),XLIM(2),XLIM(1),XLIM(1)];
                stimBoxY=[YLIM(1),YLIM(1),YLIM(2),YLIM(2),YLIM(1)];
                plot(stimBoxX,stimBoxY,'mp-','linewidth',3)
            else
%                 stimCrossX1 = [XLIM(1),XLIM(2)];
%                 stimCrossY1 = [YLIM(1),YLIM(2)];
%                 stimCrossX2 = [XLIM(1),XLIM(2)];
%                 stimCrossY2 = [YLIM(2),YLIM(1)];
                stimDotX = (XLIM(2) + XLIM(1))/2;
                stimDotY = (YLIM(2) + YLIM(1))/2;
                plot(stimDotX,stimDotY,'m.','markersize',75)
                xlim(XLIM)
                ylim(YLIM)
                set(ax,'visible','off')
                set(ax,'color','none')
            end
        end
    end
end


end