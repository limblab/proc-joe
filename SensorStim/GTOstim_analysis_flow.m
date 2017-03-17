% This script contains code to perform all GTOstim analysis in the order
% that I have developed. The goal is to be able to run a GTOstim
% experiment and then simply follow the flow of this code to perform the
% analysis. This code is NOT meant to be run all at once, but instead in
% pieces (ctrl + enter).

% Make sure to merge, sort, and split data before running this script

%% setup code -- filepath, fileprefix, generate cds files if not there, etc.
% After this section ends, all the cds files should exist :D
funcFolder = pwd;
filepath = 'D:\Lab\Data\SensorStim\Chips_20151207\';
fileprefix = 'Chips_20151207_GTOstim';
generateCDSMultipleFiles(filepath,fileprefix); % generates cds if not there
files = dir(strcat(filepath, fileprefix, '*_cds.mat'));

BIN_SIZE = 0.04;

%% for each file, pick out neurons that are "statistically" significant in 
% a bin after the stimulation starts
if((~exist(strcat(filepath,fileprefix,'_significantNeurons.mat'))))
    for i = 1:length(files)
        filename = files(i).name;
        disp(filename);
        disp(i);
        neurons{i,:} = generateAllPSTH(filepath, filename, 'noPlots',1,'useRate',1,'return',1,'binSize',BIN_SIZE);
    end
    save(strcat(filepath,fileprefix,'_significantNeurons.mat'),'neurons','files');
    disp('DONE --- can proceed');
else
    load(strcat(filepath,fileprefix,'_significantNeurons.mat'));
    disp('DONE --- can proceed');
end
%% Sometimes, it is desired to go look at neurons in files individually
% this section allows for that 
tgtFile = 2; % file number in files
filename = files(tgtFile).name;
% filename = 'Han_20170209_GTOstim_Brachialis_0175mA_area2_001_cds';
load([filepath filename]); % load in cds

neuronNumber = neurons{tgtFile}(1); % target neuron number
generatePESTH(cds, neuronNumber,'useRate',1, 'averageSpikeWaveform',12,...
    'sequenceTimes',sequenceTimes,'eventTimes',eventTimes,'confidenceInterval',1,...
    'binSize',BIN_SIZE);
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

%%%%%%%% ok analysis  code %%%%%%%%%%%%%%

%% for a neuron in common, plot all PSTH (for all muscles) -- then provide label for which muscle the neuron responds best to
neuronInCommon = 4;
plotAllMusclesForNeuron(common, neuronInCommon,'confidenceInterval',1,'label',1);
% common.neuronLabels{neuronInCommon} = input('Provide muscle label for neuron: ', 's'); % fill me in :D

%% alternatively, generate labels based on significant bins for all neurons in common
for i = 1:length(common.muscles)
    common.neuronLabels{i} = plotAllMusclesForNeuron(common, i,'confidenceInterval',1,'label',1,'noPlots',1);
end
save(strcat(filepath,fileprefix,'_significantNeuronsAllRates.mat'),'common');

%% for each neuron, plot the average relative firing rate in each 50-100ms
% bin as a scatter plot
for i = 1:length(common.muscles)
    figure();
    subplot(2,1,1)
    plotMuscleFiringRate(common,i,filepath,fileprefix,'relative',0,'figure',0);
    subplot(2,1,2)
    plotMuscleFiringRate(common,i,filepath,fileprefix,'relative',1,'figure',0,'legend',0);
end

%% attempt to label array positions based on mapfile -- need sensory mapping for this to be useful
mapfile = 'D:\Lab\Data\MapFiles\Han_Left_S1\SN 6251-001459.cmp';
arrayPositionLabel = labelArray(common,mapfile);


% for significance -- plot baseline samples per bin distribution, if normal
% then t-test, otherwise boostrap confidence interval. Also on PSTH plot,
% plot dashed line for baseline noise w/ conf interval (compute from beginning of file)

% randomely split file into two -- see if same neurons are in same file
