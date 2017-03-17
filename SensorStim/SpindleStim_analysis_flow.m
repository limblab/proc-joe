% This script contains code to perform all spindle stim analysis in the order
% that I have developed. The goal is to be able to run a spindle stim
% experiment and then simply follow the flow of this code to perform the
% analysis. This code is NOT meant to be run all at once, but instead in
% pieces (ctrl + enter).

% Make sure to merge, sort, and split data before running this script

%% setup code -- filepath, fileprefix, generate cds files if not there, etc.
% After this section ends, all the cds files should exist :D
funcFolder = pwd;
filepath = 'D:\Lab\Data\SensorStim\Han_20170217\';
fileprefix = 'Han_20170217_SpindleStim';
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
BIN_SIZE = 0.05;

%% for each file, pick out neurons that are "statistically" significant in 
% a bin after the stimulation starts

if(true || (~exist(strcat(filepath,fileprefix,'_significantNeurons.mat'))))
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
tgtFile = 9; % file number in files
filename = files(tgtFile).name;
% filename = 'Han_20170209_GTOstim_Brachialis_0175mA_area2_001_cds';
load([filepath filename]); % load in cds

neuronNumber = neurons{tgtFile}(1); % target neuron number
generatePESTH(cds, neuronNumber,'useRate',1, 'averageSpikeWaveform',0,...
    'sequenceTimes',sequenceTimes,'eventTimes',eventTimes,'confidenceInterval',0,...
    'binSize',BIN_SIZE,'plotEndStimulation',1);
% prune stuff away, use idk in neuron{tgtFile} to select
% neuronsKeep = [2];
% neurons{tgtFile} = neurons{tgtFile}(neuronsKeep);

%% Once a set of neurons have been selected for each muscle, put this info
% together in a struct called 'common'
doAgain =1; % if you want to rerun versus check and load

if(doAgain || ~exist(strcat(filepath,fileprefix,'_significantNeuronsAllRates.mat')))
    common = getAllNeuronRates(filepath, files, neurons,'binSize',BIN_SIZE);
    common.neuronLabels = common.muscles; % this is not a good labeling -- should redo in later analysis
    save(strcat(filepath,fileprefix,'_significantNeuronsAllRates.mat'),'common');
    disp('DONE -- can proceed');
else
    load(strcat(filepath,fileprefix,'_significantNeuronsAllRates.mat'));
end

%% for a neuron in common, plot all PSTH (for all muscles) -- then provide label for which muscle the neuron responds best to
neuronInCommon = 5;
plotAllMusclesForNeuron(common, neuronInCommon,'confidenceInterval',0,'label',1,'plotEndStimulation',1);
% common.neuronLabels{neuronInCommon} = input('Provide muscle label for neuron: ', 's'); % fill me in :D

%% alternatively, generate labels based on significant bins for all neurons in common
for i = 1:length(common.muscles)
    common.neuronLabels{i} = plotAllMusclesForNeuron(common, i,'confidenceInterval',1,'label',1,'noPlots',1);
end
save(strcat(filepath,fileprefix,'_significantNeuronsAllRates.mat'),'common');

