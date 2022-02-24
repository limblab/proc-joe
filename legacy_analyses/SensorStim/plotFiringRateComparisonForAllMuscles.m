function [ common ] = plotFiringRateComparisonForAllMuscles( common, neuronNumber, varargin )
% plots a comparison of the firing rate in the mean vs. the bin we care
% about for all muscles given a neuron number in common
markersize = 6;
shapes = {'o','o','s','s','^','^','+','*','x',...
    'p','p','h','h','d','d','>','<','>','<'};
filledin = [0,1,0,1,0,1,1,1,1,...
    0,1,0,1,0,1,1,1,0,0];
legendStr = common.muscleNames;

figMake = 1;
legendMake = 1;
label = 0;
labelStrings = {};
for i = 1:2:numel(varargin)
    switch varargin{i}
        case 'figMake'
            figMake = varargin{i+1};
        case 'legendMake'
            legendMake = varargin{i+1};
        case 'label'
            label = varargin{i+1};
    end
end

% Get relevant data
valsToPlot = zeros(length(common.muscleNames),5); % before stim, 40-80, after, is significant (before, after)
for i = 1:length(common.muscleNames) % for each muscle
    muscleData = common.(common.muscleNames{i});
    zeroIdx=find(muscleData.edges{neuronNumber} == 0);
%     valsToPlot(i,1) = mean([muscleData.counts{neuronNumber}(1:zeroIdx),muscleData.counts{neuronNumber}(zeroIdx+2:end)]);
    valsToPlot(i,1) = mean(muscleData.counts{neuronNumber}(1:zeroIdx));
    valsToPlot(i,2) = muscleData.counts{neuronNumber}(zeroIdx+1);
    valsToPlot(i,3) = mean(muscleData.counts{neuronNumber}(zeroIdx+2:end));   
    
    % paired t-test to see if significantly different than both. mean data =
    valsToPlot(i,4) = 0;
    valsToPlot(i,5) = 0;
    numStims = numel(muscleData.eventTimes);
    t(1,1) = (valsToPlot(i,2)-valsToPlot(i,1))/(muscleData.stdDiff(neuronNumber,1)/sqrt(numStims));
    t(1,2) = (valsToPlot(i,2)-valsToPlot(i,3))/(muscleData.stdDiff(neuronNumber,2)/sqrt(numStims));
    if(1-tcdf(t(1,1),numStims-1) < 0.05/length(common.muscleNames))
        valsToPlot(i,4) = 1;
    end
    if(1-tcdf(t(1,1),numStims-1) < 0.05/length(common.muscleNames))
        valsToPlot(i,5) = 1;
    end
        
end
% plot data
if(figMake)
    figure();
else
    
hold on
greyScaleColors = linspace(0,0.35,length(common.muscleNames));
for i = 1:length(common.muscleNames)    
    % check if "significantly" better
    muscleData = common.(common.muscleNames{i});
    binCounts = muscleData.counts{neuronNumber};
    binEdges = muscleData.edges{neuronNumber};
    binSize = binEdges(2)-binEdges(1);
    zeroIdx=find(muscleData.edges{neuronNumber} == 0);
    plusMinus = muscleData.bootstrapBounds(neuronNumber,:);
    
    if(binCounts(zeroIdx+1) > plusMinus(2)/binSize && ... % bootstrap
            binCounts(zeroIdx+1)>mean(binCounts)+std(binCounts))
%     if(valsToPlot(i,4) == 1 && valsToPlot(i,5) == 1)
        % add to labelStrings, plot later
        labelStrings{end+1} = common.muscleNames{i};
    else
        colorPH = repmat(greyScaleColors(i),1,3);
        ph = plot([0 1 2], valsToPlot(i,1:3),'linewidth',2.25,'marker',shapes{i},'markersize',markersize,...
            'markerfacecolor',colorPH,'markeredgecolor',colorPH,'color',colorPH);  
        if(~filledin(i))
            set(ph,'markerfacecolor','none');
        end
    end
end

for i = 1:length(common.muscleNames)    
    % check if "significantly" better
    muscleData = common.(common.muscleNames{i});
    binCounts = muscleData.counts{neuronNumber};
    binEdges = muscleData.edges{neuronNumber};
    binSize = binEdges(2)-binEdges(1);
    zeroIdx=find(muscleData.edges{neuronNumber} == 0);
    plusMinus = muscleData.bootstrapBounds(neuronNumber,:);
%     
    if(binCounts(zeroIdx+1) > plusMinus(2)/binSize && ... % bootstrap
            binCounts(zeroIdx+1)>mean(binCounts)+std(binCounts))
       
%     if(valsToPlot(i,4) == 1 && valsToPlot(i,5) == 1) 
        % plot these now so that they are on top of the order
        c='r';
        ph = plot([0 1 2], valsToPlot(i,1:3),'linewidth',2.25,'marker',shapes{i},'markersize',markersize,...
            'markerfacecolor',c,'markeredgecolor',c,'color',c);
        if(~filledin(i))
            set(ph,'markerfacecolor','none');
        end
    end
    
end
a = gca;
a.YTick = a.YLim;
a.YLim(2) = a.YLim(2) + 1;
a.XLim = [-0.1 2.1];
a.XTick = [0, 1, 2];
a.XTickLabel = {'Before Stim','40-80ms','After Stim'};
ylabel('Average Firing Rate (spikes/s)')
formatForLee(gca);
a.XMinorTick = 'off';

if(legendMake)
    for i = 1:length(common.muscleNames)
        h(i)=scatter(-100,-100,'marker',shapes{i},'markeredgecolor','k');
        if(filledin(i))
            set(h(i),'markerfacecolor','k');
        end
    end
    a.XLim = [-0.1 2.1];
    leg = legend(h,legendStr);
end

if(label)
    if(iscell(labelStrings))
        common.neuronLabels{neuronNumber,1} = labelStrings;
    else
        common.neuronLabels{neuronNumber,1} = {labelStrings};
    end
end


end

