function [ neuronLabel ] = plotAllMusclesForNeuron(common, neuronInCommon,varargin)

confInt = 0;
label = 0;
noPlots =0;
useRate = 1;
plotEndStimulation = 1;
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
    figure();
end
    
for i = 1:length(common.muscleNames)
    bE = common.(genvarname(common.muscleNames{i})).edges{neuronInCommon};
    diffEdge = bE(2)-bE(1);
    bEplot = bE(1:end-1)+diffEdge/2;
    bC = common.(genvarname(common.muscleNames{i})).counts{neuronInCommon};
    if(~noPlots)
        subplot(numRows,numCols,i);
        bar(bEplot,bC);
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

            h1=plot([-100,100],[meanConf,meanConf],'g','linewidth',2)
            plot([-100,100],[plusMinus(1),plusMinus(1)],'g','linewidth',2)
            plot([-100,100],[plusMinus(2),plusMinus(2)],'g','linewidth',2)

            meanCounts = mean(bC);
            stdCounts = std(bC);
            plot([-100,100],[meanCounts,meanCounts],'r','linewidth',2);
            h2=plot([-100,100],[meanCounts+2*stdCounts,meanCounts+2*stdCounts],'r','linewidth',2);

            if(i==1)
                legend([h1 h2],'Boostrap','Mean+2*Std');
            end
            a.XLim = xLimits;
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

