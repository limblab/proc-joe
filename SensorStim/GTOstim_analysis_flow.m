% This script contains code to perform all GTOstim analysis in the order
% that I have developed. The goal is to be able to run a GTOstim
% experiment and then simply follow the flow of this code to perform the
% analysis. This code is NOT meant to be run all at once, but instead in
% pieces (ctrl + enter).

% Make sure to merge, sort, and split data before running this script

%% setup code -- filepath, fileprefix, generate cds files if not there, etc.
% After this section ends, all the cds files should exist :D
funcFolder = pwd;
filepath = 'D:\Lab\Data\SensorStim\Han_20170209\';
fileprefix = 'Han_20170209_GTOstim';
generateCDSMultipleFiles(filepath,fileprefix); % generates cds if not there
files = dir(strcat(filepath, fileprefix, '*_cds.mat'));
% Han_20170209
% Han_20170217

BIN_SIZE = 0.04;

%% for each file, pick out neurons that are "statistically" significant in 
% a bin after the stimulation starts
if((~exist(strcat(filepath,fileprefix,'_significantNeurons.mat'))))
    for i = 1:length(files)
        filename = files(i).name;
        disp(filename);
        disp(i);
        neurons{i,:} = generateAllPSTH(filepath, filename, 'noPlots',1,'useRate',1,'return',1,'binSize',BIN_SIZE,...
            'spindleStim',0);
    end
    save(strcat(filepath,fileprefix,'_significantNeurons.mat'),'neurons','files');
    disp('DONE --- can proceed');
else
    load(strcat(filepath,fileprefix,'_significantNeurons.mat'));
    disp('DONE --- can proceed');
end
%% Sometimes, it is desired to go look at neurons in files individually
% this section allows for that 
figDir = 'D:\Lab\Summary Figures\SensorStim\GTOstim\';
for tgtFile = 1:size(neurons,1)% file number
    filename = files(tgtFile).name;
    idxUnderscore = strfind(filename,'_');
    muscleName = filename(idxUnderscore(3)+1:idxUnderscore(4)-1);
    dateRan = filename(idxUnderscore(1)+1:idxUnderscore(2)-1);
    load([filepath filename]); % load in cds
    figName = strcat(filename(1:end-8),'_');
    figDirCurr = strcat(figDir,filesep,muscleName,'_',dateRan);
    for nn = 1:numel(neurons{tgtFile,1})
        neuronNumber = neurons{tgtFile}(nn); % target neuron number
        generatePESTH(cds, neuronNumber,'useRate',1, 'averageSpikeWaveform',12,...
            'sequenceTimes',sequenceTimes,'eventTimes',eventTimes,'confidenceInterval',1,...
            'binSize',BIN_SIZE,'binsAbove',0,'highlightBin',1);
        
        f = figure(1);
        formatForLee(f);
        saveFigure(f,figDirCurr,strcat(figName,'nn',num2str(neuronNumber),'_histogram'));
        f = figure(2);
        formatForLee(f);
        saveFigure(f,figDirCurr,strcat(figName,'nn',num2str(neuronNumber),'_averageWaveform'));
        close all
    end
    disp(tgtFile)
end
% 
for tgtFile = 1:size(neurons,1)% file number
    filename = files(tgtFile).name;
    idxUnderscore = strfind(filename,'_');
    muscleName = filename(idxUnderscore(3)+1:idxUnderscore(4)-1);
    dateRan = filename(idxUnderscore(1)+1:idxUnderscore(2)-1);
    load([filepath filename]); % load in cds
    figName = strcat(filename(1:end-8),'_');
    figDirCurr = strcat(figDir,filesep,muscleName,'_',dateRan);
    for nn = 1:numel(neurons{tgtFile,1})
        neuronNumber = neurons{tgtFile}(nn); % target neuron number
        generateRaster(cds,neuronNumber,eventTimes,sequenceTimes,'plotStimTime',1,'GTOstim',1);
        formatForLee(gcf);
        saveFigure(gcf,figDirCurr,strcat(figName,'nn',num2str(neuronNumber),'_raster'));
        close all
    end
    disp(tgtFile)
end

%% prune stuff away, use idk in neuron{tgtFile} to select
neuronsKeep = [2];
neurons{tgtFile} = neurons{tgtFile}(neuronsKeep);
%% Once a set of neurons have been selected for each muscle, put this info
% together in a struct called 'common'
doAgain = 1; % if you want to rerun versus check and load

if(doAgain || ~exist(strcat(filepath,fileprefix,'_significantNeuronsAllRates.mat')))
    common = getAllNeuronRates(filepath, files, neurons,'binSize',BIN_SIZE);
    common.neuronLabels = common.muscles; % this is not a good labeling -- should redo in later analysis
    common.fileprefix = fileprefix;
    common.filepath = filepath;
    save(strcat(filepath,fileprefix,'_significantNeuronsAllRates.mat'),'common');
    disp('DONE -- can proceed');
else
    load(strcat(filepath,fileprefix,'_significantNeuronsAllRates.mat'));
end

%%%%%%%% ok analysis  code %%%%%%%%%%%%%%

%% for a neuron in common, plot all PSTH (for all muscles) -- then provide label for which muscle the neuron responds best to
neuronInCommon = 3;
plotAllMusclesForNeuron(common, neuronInCommon,'confidenceInterval',1,'label',1);
% common.neuronLabels{neuronInCommon} = input('Provide muscle label for neuron: ', 's'); % fill me in :D

%% for each neuron or a neuron, plot change in firing rate in bin vs. mean for each muscle
figHandle = figure;
% [p,] = numSubplots(length(common.muscles));
for i = 1:length(common.muscles)
%     subplot(p(1),p(2),i);
    if(i==1)
        legendMake = 1;
    else
        legendMake = 0;
    end
    common=plotFiringRateComparisonForAllMuscles(common,i,'figMake',0,'legendMake',legendMake,'label',1);
    save(strcat(filepath,fileprefix,'_significantNeuronsAllRates.mat'),'common'); 
    close all
    %     title(strcat('NN: ', num2str(i)))
end

%% for each muscle, plot change in firing rate in bin vs. mean for all neurons
% figHandle = figure;
% [p,] = numSubplots(length(common.muscleNames));
% for i = 1:length(common.muscleNames)
%     subplot(p(1),p(2),i);
%     if(i==1)
%         legendMake = 0;
%     else
%         legendMake = 0;
%     end
%     plotFiringRateComparisonForAllNeurons(common,common.muscleNames{i},'figMake',0,'legendMake',legendMake);
%     title(common.muscleNames{i})
% end

%% bar plots of significant neuron per each muscle
musclesPerNeuron = categorical(common.muscleNames);
neuronPerMuscleCount = zeros(1,numel(common.muscleNames));
for i = 1:numel(common.neuronLabels)
    for j = 1:numel(common.neuronLabels{i})
        for mn = 1:numel(common.muscleNames)
            if(iscell(common.neuronLabels{i}))
                if(strcmp(common.neuronLabels{i}(j),common.muscleNames{mn})==1)
                    neuronPerMuscleCount(1,mn) = neuronPerMuscleCount(1,mn)+1;
                end
            else
                if(strcmp(common.neuronLabels{i},common.muscleNames{mn})==1)
                    neuronPerMuscleCount(1,mn) = neuronPerMuscleCount(1,mn)+1;
                end
            end
            
        end
    end
end
b=bar(neuronPerMuscleCount);
ylabel('# Neurons driven')
set(gca,'xtick',[1:numel(common.muscleNames)])
set(gca,'xticklabel',common.muscleNames)
set(gca,'XTickLabelRotation',45)
formatForLee(gcf);
set(gca, 'XMinorTick', 'off')
set(gca,'XTick',[1:numel(common.muscleNames)]);
set(gca,'YMinorTick','off');
set(gca,'YTick',0:max(neuronPerMuscleCount));

%% bar plots of significant muscle number per each neuron
musclesPerNeuronCount = zeros(1,numel(common.neuronLabels));
for i = 1:numel(common.neuronLabels)
    musclesPerNeuronCount(i) = numel(common.neuronLabels{i});
end
b=bar(musclesPerNeuronCount);
xlabel('Neuron Number')
ylabel('# muscles')
formatForLee(gcf);
set(gca,'Xtick',[1:numel(common.muscleNames)])
set(gca, 'XMinorTick', 'off')
set(gca,'YMinorTick','off');
set(gca,'YTick',0:max(neuronPerMuscleCount));
xlim([0.2,4.8])


%% Compare across multiple files
filePrefixesToCompare = {'Chips_20151208_GTOstim'};
filePathsToCompare = {'D:\Lab\Data\SensorStim\Chips_20151208\'};
compareAcrossFiles(common,2,filePathsToCompare,filePrefixesToCompare);
disp('Plotted!');
%% for each neuron, plot the average relative firing rate in each 50-100ms
% bin as a scatter plot
% for i = 1:length(common.muscles)
%     figure();
%     subplot(2,1,1)
%     plotMuscleFiringRate(common,i,filepath,fileprefix,'relative',0,'figure',0);
%     subplot(2,1,2)
%     plotMuscleFiringRate(common,i,filepath,fileprefix,'relative',1,'figure',0,'legend',0);
% end

% for significance -- plot baseline samples per bin distribution, if normal
% then t-test, otherwise boostrap confidence interval. Also on PSTH plot,
% plot dashed line for baseline noise w/ conf interval (compute from beginning of file)

% randomely split file into two -- see if same neurons are in same file
