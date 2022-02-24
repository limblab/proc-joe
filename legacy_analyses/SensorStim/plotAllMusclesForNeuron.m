function [ neuronLabel ] = plotAllMusclesForNeuron(common, neuronInCommon,varargin)

confInt = 0;
label = 0;
noPlots =0;
useRate = 1;
plotEndStimulation = 0;
legendFontSize = 14;
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'confidenceInterval'
            confInt = varargin{i+1};
        case 'label'
            label = varargin{i+1};
        case 'noPlots'
            noPlots = varargin{i+1};
        case 'useRate'
            useRate = varargin{i+1};
        case 'plotEndStimulation'
            plotEndStimulation = varargin{i+1};

    end
end

if(label)
    neuronLabel = {};
end

numPlots = length(common.muscleNames);
[p,] = numSubplots(numPlots);
numCols = p(2);
numRows = p(1);

if(~noPlots)
    barFig = figure();
end
    
for i = 1:length(common.muscleNames)
    bE = common.(genvarname(common.muscleNames{i})).edges{neuronInCommon};
    binSize = bE(2)-bE(1);
    bEplot = bE(1:end-1)+binSize/2;
    bC = common.(genvarname(common.muscleNames{i})).counts{neuronInCommon};
    if(~noPlots)
        subplot(numRows,numCols,i);
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
        set(gca,'XLim',aHand.XLim,'YLim',aHand.YLim)
        xlabel('Time After Stimulation (s)')
        ylabel('Average Firing Rate (spikes/s)');
        if(confInt)
            a =  gca;
            xLimits = a.XLim;
            hold on
            meanConf = common.(genvarname(common.muscleNames{i})).bootstrapMean(neuronInCommon,1);
            plusMinus = common.(genvarname(common.muscleNames{i})).bootstrapBounds(neuronInCommon,:);
            if(useRate)
                meanConf = meanConf./(bE(2)-bE(1));
                plusMinus = plusMinus./(bE(2)-bE(1));
            end
            
            h1=plot([binSize,2*binSize],[plusMinus(2),plusMinus(2)],'g','linewidth',2);

            meanCounts = mean(bC);
            stdCounts = std(bC);
            h2=plot([binSize,2*binSize],[meanCounts+2*stdCounts,meanCounts+2*stdCounts],'b','linewidth',2);
            
            binSize = bE(2)-bE(1);
            eventTimes = common.(common.muscleNames{i}).eventTimes;
%             if(useRate)
%                 poissVal = poissinv(0.975,meanCounts*(binSize)*numel(eventTimes))/numel(eventTimes)/binSize;
%             else
%                 poissVal = poissinv(0.975,meanCounts*numel(eventTimes))/numel(eventTimes);
%             end
%             h3 = plot([-100,100],[poissVal, poissVal],'m','linewidth',2);
%             
            if(i==1)
                l=legend([h1 h2],'Boostrap','Mean+2*Std');
                set(l,'box','off','fontsize',legendFontSize)
            end
            a.XLim = xLimits;
            formatForLee(a);
        end      
        if(plotEndStimulation)
            stimSampRate = 1/1e-3;
            
            stimStateTimeAverage = (sum(common.(genvarname(common.muscleNames{i})).stimState == 1))/numel(common.(genvarname(common.muscleNames{i})).eventTimes)/stimSampRate;

            % plot vertical line at this time
            hold on
            figHandle = gcf;
            plot([stimStateTimeAverage, stimStateTimeAverage],figHandle.CurrentAxes.YLim,'--k','linewidth',2)
        end
        title(common.muscleNames{i});
    end
    
    if(label && keepNeuronGTOstim(-1,-1,bE,bC,common.(genvarname(common.muscleNames{i})).eventTimes(1),1,...
            'bootstrap',0,'meanCounts',common.(genvarname(common.muscleNames{i})).bootstrapMean(neuronInCommon,1),...
            'plusMinus',common.(genvarname(common.muscleNames{i})).bootstrapBounds(neuronInCommon,:))==1)  
        
        neuronLabel{end+1,1} = common.muscleNames{i};
    end
end


end

