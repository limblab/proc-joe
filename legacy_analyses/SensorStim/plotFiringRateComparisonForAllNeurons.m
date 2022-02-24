function [  ] = plotFiringRateComparisonForAllNeurons( common, muscleName, varargin )
% plots a comparison of the firing rate in the mean vs. the bin we care
% about for all muscles given a neuron number in common
markersize = 6;

shapes = {'o','o','s','s','^','^','+','*','x',...
    'p','p','h','h','d','d','>','<','>','<','.','.'};
filledin = [0,1,0,1,0,1,1,1,1,...
    0,1,0,1,0,1,1,1,0,0,0,0];
legendStr = {'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21'};
figMake = 1;
legendMake = 1;
for i = 1:2:numel(varargin)
    switch varargin{i}
        case 'figMake'
            figMake = varargin{i+1};
        case 'legendMake'
            legendMake = varargin{i+1};
    end
end

% Get relevant data
valsToPlot = zeros(length(common.muscles),2);
for i = 1:length(common.muscles) % for each neuron
    muscleData = common.(muscleName);
    zeroIdx=find(muscleData.edges{i} == 0);
    valsToPlot(i,1) = mean([muscleData.counts{i}(1:zeroIdx),muscleData.counts{i}(zeroIdx+2:end)]);
%     valsToPlot(i,1) = mean(muscleData.counts{neuronNumber}(1:zeroIdx));
    valsToPlot(i,2) = muscleData.counts{i}(zeroIdx+1);
%     valsToPlot(i,3) = mean(muscleData.counts{neuronNumber}(zeroIdx+2:end));   
end
slopes = (valsToPlot(:,2) - valsToPlot(:,1));
% plot data 
if(figMake)
    figure();
else
    
hold on
greyScaleColors = linspace(0,0.35,length(common.muscles));
for i = 1:length(common.muscles)    
    % check if "significantly" better
    muscleData = common.(muscleName);
    binCounts = muscleData.counts{i};
    binEdges = muscleData.edges{i};
    binSize = binEdges(2)-binEdges(1);
    zeroIdx=find(muscleData.edges{i} == 0);
    plusMinus = muscleData.bootstrapBounds(i,:);
    poissVal = poissinv(0.975,mean(binCounts)*(binSize)*numel(muscleData.eventTimes))/numel(muscleData.eventTimes)/binSize;
    
    if(binCounts(zeroIdx+1) > plusMinus(2))
        % do nothing yet, plot all of these at the end
    else
        colorPH = repmat(greyScaleColors(i),1,3);
        ph = plot([0 1], valsToPlot(i,:),'linewidth',2.25,'marker',shapes{i},'markersize',markersize,...
            'markerfacecolor',colorPH,'markeredgecolor',colorPH,'color',colorPH);  
        if(~filledin(i))
            set(ph,'markerfacecolor','none');
        end
    end
end

for i = 1:length(common.muscles)    
    % check if "significantly" better
    muscleData = common.(muscleName);
    binCounts = muscleData.counts{i};
    binEdges = muscleData.edges{i};
    binSize = binEdges(2)-binEdges(1);
    zeroIdx=find(muscleData.edges{i} == 0);
    plusMinus = muscleData.bootstrapBounds(i,:);
    poissVal = poissinv(0.975,mean(binCounts)*(binSize)*numel(muscleData.eventTimes))/numel(muscleData.eventTimes)/binSize;
    
    if(binCounts(zeroIdx+1) > plusMinus(2))
        % plot these now so that they are on top of the order
        c='r';
        ph = plot([0 1], valsToPlot(i,:),'linewidth',2.25,'marker',shapes{i},'markersize',markersize,...
            'markerfacecolor',c,'markeredgecolor',c,'color',c);
        if(~filledin(i))
            set(ph,'markerfacecolor','none');
        end
    end
    
end
a = gca;
a.YTick = a.YLim;
a.YLim(2) = a.YLim(2) + 1;
a.XLim = [-0.1 1.1];
a.XTick = [0, 1];
a.XTickLabel = {'Average','After Stim'};
formatForLee(gca);
a.XMinorTick = 'off';

if(legendMake)
    for i = 1:length(common.muscles)
        h(i)=scatter(-100,-100,'marker',shapes{i},'markeredgecolor','k');
        if(filledin(i))
            set(h(i),'markerfacecolor','k');
        end
    end
    a.XLim = [-0.1 1.1];
    leg = legend(h,legendStr{1:length(common.muscles)});
end


end

