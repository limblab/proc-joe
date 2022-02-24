function [ figHandle ] = plotMuscleFiringRate( common, num,filepath,fileprefix, varargin )

markersize = 15;
% shape, filled in, muscle name dict
shapes = {'o','o','s','o','o','s','o','o','s','o','o','s','o','o','s','o','o','s','o','o','s'};
filledin = {0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1};
colors = {'b','b','b','r','r','r','g','g','g','k','k','k','m','m','m','c','c','c','y','y','y'};
muscleName = {'ADELT','BRAD','Brachialis','ECR','ECU','EDC','EDChigh','EDClow','FCR',...
    'FCU','FDS','InfPec','PDELT2','PDELT','SupPec','TRIlat','TRImed','medBI'};
% plots all average firing rates (relative if requested) for the given neuron across all muscles

relative = 0;
bin = 1;
figMake = 1;
leg = 1;
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'relative'
            relative = varargin{i+1};
        case 'bin'
            bin = varargin{i+1};
        case 'figure'
            figMake = varargin{i+1};
        case 'legend'
            leg = varargin{i+1};
    end
end

% get all vals to plot
vals = zeros(length(common.muscleNames),1);
for i = 1:length(common.muscleNames)
    zeroIdx = find(common.(genvarname(common.muscleNames{i})).edges{num}==0);
    idxBin = zeroIdx + bin;
    c = common.(genvarname(common.muscleNames{i})).counts{num};
    vals(i,1) = c(1,idxBin);
    if(relative)
%         % load appropriate cds
%         f = dir(strcat(filepath,fileprefix,'_*',common.muscleNames{i},'*_cds.mat'));
%         load(strcat(filepath,f(1).name));
%         % do bootstrap
%         neuronNumber = findNeuronNumber(cds,common.channels(num),common.IDs(num),common.electrodes{num});
%         [divisor, plusMinus] = bootstrapConfidenceInterval(cds, neuronNumber, eventTimes(1), 0.05);
%         
        % divide
        divisor = mean(c(1,:));
        vals(i,1) = vals(i,1)/divisor;
        
        % clear cds
        clear cds; clear eventTimes; clear sequenceTimes;
    end
end

if(figMake)
    figure();
end

set(gca,'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
         'YLim',[0 eps],...             %#   set the y axis limit (tiny!),
         'Color','none');
hold on;   
% hAxes = axes('NextPlot','add',...           %# Add subsequent plots to the axes,
%          'DataAspectRatio',[1 1 1],...  %#   match the scaling of each axis,
%          'YLim',[0 eps],...             %#   set the y axis limit (tiny!),
%          'Color','none');               %#   and don't use a background color

for i = 1:length(vals)    
    if(filledin{i} == 1)
        plot(vals(i,1),0,'marker',shapes{i},'MarkerFaceColor',colors{i},'MarkerEdgeColor',...
            colors{i},'markersize',markersize-i*0.5)
    else
        plot(vals(i,1),0,'marker',shapes{i},'MarkerFaceColor','w','MarkerEdgeColor',...
            colors{i},'markersize',markersize-i*0.5)
    end
end
if(leg)
    legend(muscleName);
end
end

