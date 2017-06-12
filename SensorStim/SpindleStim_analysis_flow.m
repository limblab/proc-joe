% This script contains code to perform all spindle stim analysis in the order
% that I have developed. The goal is to be able to run a spindle stim
% experiment and then simply follow the flow of this code to perform the
% analysis. This code is NOT meant to be run all at once, but instead in
% pieces (ctrl + enter).

% Make sure to merge, sort, and split data before running this script

%% setup code -- filepath, fileprefix, generate cds files if not there, etc.
% After this section ends, all the cds files should exist :D
funcFolder = pwd;
filepath = 'D:\Lab\Data\SensorStim\Lando_20170511\';
fileprefix = 'Lando_20170511_SpindleStim';
generateCDSMultipleFiles(filepath,fileprefix); % generates cds if not there
files = dir(strcat(filepath, fileprefix, '*_cds.mat'));
% i = 1;
% while i <= length(files)
%     if(~isempty(strfind(files(i).name,'-s')))
%         files(i) = []; % removes row i from struct files
%         i=i-1;
%     elseif(~isempty(strfind(files(i).name,'GTO')))
%         files(i) = []; % removes row i from struct files
%         i=i-1;
%     elseif(~isempty(strfind(files(i).name,'Sweep')))
%         files(i) = []; % removes row i from struct files
%         i=i-1;
%     elseif(~isempty(strfind(lower(files(i).name),'neurons')))
%         files(i) = [];
%         i=i-1;
%     end
%     i=i+1;
% end
BIN_SIZE = 0.04;

%% for each file, pick out neurons that are "statistically" significant in 
% a bin after the stimulation starts

if((~exist(strcat(filepath,fileprefix,'_significantNeurons.mat'))))
    for i = 1:length(files)
        filename = files(i).name;
        disp(filename);
        disp(i);
        neurons{i,:} = generateAllPSTH(filepath, filename, 'noPlots',1,'useRate',1,'return',1,'binSize',BIN_SIZE,...
            'spindleStim',1);
    end
    save(strcat(filepath,fileprefix,'_significantNeurons.mat'),'neurons','files');
    disp('DONE --- can proceed');
else
    load(strcat(filepath,fileprefix,'_significantNeurons.mat'));
    disp('DONE --- can proceed');
end

%% Sometimes, it is desired to go look at neurons in files individually
% this section allows for that 

figDir = 'D:\Lab\Summary Figures\SensorStim\Spindlestim\';
for tgtFile = 1%:size(neurons,1)% file number
    filename = files(tgtFile).name;
    idxUnderscore = strfind(filename,'_');
    muscleName = filename(idxUnderscore(3)+1:idxUnderscore(4)-1);
    dateRan = filename(idxUnderscore(1)+1:idxUnderscore(2)-1);
    load([filepath filename]); % load in cds
    figName = strcat(filename(1:end-8),'_');
    figDirCurr = strcat(figDir,filesep,muscleName,'_',dateRan);
    for nn  = 1:numel(neurons{tgtFile,1})
        neuronNumber = neurons{tgtFile}(nn); % target neuron number
        [bE,bC,~]=generatePESTH(cds, neuronNumber,'spindleStim',1,'useRate',1, 'averageSpikeWaveform',0,...
            'sequenceTimes',sequenceTimes,'eventTimes',eventTimes,'stimState',stimState,'confidenceInterval',1,...
            'binSize',0.01,'plotEndStimulation',1,'highlightBin',0);
        formatForLee(gcf)
%         saveFigure(gcf,figDirCurr,strcat(figName,'nn',num2str(neuronNumber),'_histogram'));
%         close all
    end
    disp(tgtFile)
end
% 
% for tgtFile = 1:size(neurons,1)% file number
%     filename = files(tgtFile).name;
%     idxUnderscore = strfind(filename,'_');
%     muscleName = filename(idxUnderscore(3)+1:idxUnderscore(4)-1);
%     dateRan = filename(idxUnderscore(1)+1:idxUnderscore(2)-1);
%     load([filepath filename]); % load in cds
%     figName = strcat(filename(1:end-8),'_');
%     figDirCurr = strcat(figDir,filesep,muscleName,'_',dateRan);
%     for nn = 1:numel(neurons{tgtFile,1})
%         neuronNumber = neurons{tgtFile}(nn); % target neuron number
%         generateRaster(cds,neuronNumber,eventTimes,sequenceTimes,'plotStimTime',1,'GTOstim',0);
% %         saveFigure(gcf,figDirCurr,strcat(figName,'nn',num2str(neuronNumber),'_raster'));
% %         close all
%     end
%     disp(tgtFile)
% end

%% Once a set of neurons have been selected for each muscle, put this info
% together in a struct called 'common'
doAgain =0; % if you want to rerun versus check and load

if(doAgain || ~exist(strcat(filepath,fileprefix,'_significantNeuronsAllRates.mat')))
    common = getAllNeuronRates(filepath, files, neurons,'binSize',BIN_SIZE);
    common.neuronLabels = common.muscles; % this is not a good labeling -- should redo in later analysis
    save(strcat(filepath,fileprefix,'_significantNeuronsAllRates.mat'),'common');
    disp('DONE -- can proceed');
else
    load(strcat(filepath,fileprefix,'_significantNeuronsAllRates.mat'));
end

%% for a neuron in common, plot all PSTH (for all muscles) -- then provide label for which muscle the neuron responds best to
neuronInCommon = 1;
plotAllMusclesForNeuron(common, neuronInCommon,'confidenceInterval',0,'label',1,'plotEndStimulation',1);
% common.neuronLabels{neuronInCommon} = input('Provide muscle label for neuron: ', 's'); % fill me in :D

%% alternatively, generate labels based on significant bins for all neurons in common
common.neuronLabels = {};
for i = 1:length(common.muscles)
    common.neuronLabels{i} = {};
    common.neuronLabels{i} = plotAllMusclesForNeuron(common, i,'noPlots',1,'label',1);
end
save(strcat(filepath,fileprefix,'_significantNeuronsAllRates.mat'),'common');

%% Look at EMG data (if exists) 
for i = 2:2%numel(files)
    % load cds
    load(strcat(filepath,files(i).name));
    % call function
    spindleStimEMGAnalysis(cds,eventTimes,files(i).name);
end