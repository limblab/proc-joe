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

BIN_SIZE = 0.05;

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
tgtFile = 9; % file number in files
filename = files(tgtFile).name;
% filename = 'Han_20170209_GTOstim_Brachialis_0175mA_area2_001_cds';
load([filepath filename]); % load in cds

neuronNumber = neurons{tgtFile}(1); % target neuron number
generatePESTH(cds, neuronNumber,'useRate',1, 'averageSpikeWaveform',0,...
    'sequenceTimes',sequenceTimes,'eventTimes',eventTimes,'confidenceInterval',1,...
    'binSize',BIN_SIZE,'plotEndStimulation',1);
% prune stuff away, use idk in neuron{tgtFile} to select
% neuronsKeep = [2];
% neurons{tgtFile} = neurons{tgtFile}(neuronsKeep);