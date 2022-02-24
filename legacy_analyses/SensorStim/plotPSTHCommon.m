function [  ] = plotPSTHCommon( common, neuronInCommon,muscleName, varargin )

confInt = 0;
noPlots =0;
useRate = 1;
plotEndStimulation = 0;
figMake = 1;
legendMake = 1;
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'confidenceInterval'
            confInt = varargin{i+1};
        case 'noPlots'
            noPlots = varargin{i+1};
        case 'useRate'
            useRate = varargin{i+1};
        case 'plotEndStimulation'
            plotEndStimulation = varargin{i+1};
        case 'figMake'
            figMake = varargin{i+1};
        case 'legendMake'
            legendMake = varargin{i+1};
    end
end

if(~noPlots && figMake)
    barFig = figure();
else
    barFig = gcf;
end

bE = common.(muscleName).edges{neuronInCommon};
binSize = bE(2)-bE(1);
bEplot = bE(1:end-1)+binSize/2;
bC = common.(muscleName).counts{neuronInCommon};
if(~noPlots)
    barWidth = 0.95;
    bf = bar(bEplot,bC,'edgecolor','none','facecolor','k','barWidth',barWidth);
    % color bar we care about a different color
    aHand = barFig.CurrentAxes;
    zeroIdx = find(bE==0);
    xData = aHand.Children.XData;
    yData = aHand.Children.YData*0;
    yData(zeroIdx+1) = aHand.Children.YData(zeroIdx+1);
    hold on
    bf = bar(xData,yData,...
        'parent',aHand,'faceColor','r','EdgeColor','none','barWidth',barWidth);
    set(gca,'XLim',aHand.XLim,'YLim',aHand.YLim);
    xlabel('Time Relative to Start (s)');
    ylabel('Average Firing Rate (spikes/s)');
    
    % color bar after zero a different color
    aHand = barFig.CurrentAxes;
    zeroIdx = find(bE==0);
    xData = aHand.Children(2).XData;
    yData = aHand.Children(2).YData*0;
    yData(zeroIdx) = aHand.Children(2).YData(zeroIdx);
    hold on
    bf = bar(xData,yData,...
        'parent',aHand,'faceColor',[0.25 0.25 0.25],'EdgeColor','none','barWidth',barWidth);
    set(gca,'XLim',aHand.XLim,'YLim',aHand.YLim);
    xlabel('Time Relative to Start (s)');
    ylabel('Average Firing Rate (spikes/s)');
    
    
    if(confInt)
        a =  gca;
        xLimits = a.XLim;
        hold on
        meanConf = common.(muscleName).bootstrapMean(neuronInCommon,1);
        plusMinus = common.(muscleName).bootstrapBounds(neuronInCommon,:);
        if(useRate)
            meanConf = meanConf./(bE(2)-bE(1));
            plusMinus = plusMinus./(bE(2)-bE(1));
        end

        h1=plot([-100,2*100],[plusMinus(2),plusMinus(2)],'g','linewidth',2);

        meanCounts = mean(bC);
        stdCounts = std(bC);
        h2=plot([-100,100],[meanCounts+2*stdCounts,meanCounts+2*stdCounts],'b','linewidth',2);

        binSize = bE(2)-bE(1);
        eventTimes = common.(muscleName).eventTimes;
%         if(useRate)
%             poissVal = poissinv(0.975,meanCounts*(binSize)*numel(eventTimes))/numel(eventTimes)/binSize;
%         else
%             poissVal = poissinv(0.975,meanCounts*numel(eventTimes))/numel(eventTimes);
%         end
%         h3 = plot([-100,100],[poissVal, poissVal],'m','linewidth',2);
%           
        if(legendMake)
            l=legend([h1 h2],'Boostrap','Mean+2*Std');
            set(l,'box','off','fontsize',12)
        end
        a.XLim = xLimits;
        a.YLim(2) = max(a.YLim(2),plusMinus(2))*1.3;
    end      
    if(plotEndStimulation)
        stimSampRate = 1/1e-3;

        stimStateTimeAverage = (sum(common.(genvarname(common.muscleNames{i})).stimState == 1))/numel(common.(genvarname(common.muscleNames{i})).eventTimes)/stimSampRate;

        % plot vertical line at this time
        hold on
        figHandle = gcf;
        plot([stimStateTimeAverage, stimStateTimeAverage],figHandle.CurrentAxes.YLim,'--k','linewidth',2)
    end
end
formatForLee(barFig);

end

