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
fileprefix = 'Han_20170209';
generateCDSMultipleFiles(filepath,fileprefix); % generates cds if not there

%% for each file, pick out neurons that are "statistically" significant in 
% a bin after the stimulation starts

files = dir(strcat(filepath, fileprefix, '*_cds.mat'));

for i = 1:length(files)
    filename = files(i).name;
    disp(filename);
    disp(i);
    generateAllPSTH(filepath, filename);
    pause; % to give time for user to do any documenting
end

%% Sometimes, it is desired to go look at neurons in files individually
% this section allows for that 
% tgtFile = 4; % file number in files
% filename = files(tgtFile).name;
filename = 'Han_20170209_GTOstim_ECU_08mA_area2_001_cds';
load([filepath filename]); % load in cds

neuronNumber =73; % target neuron number
generatePESTH(cds, neuronNumber,'useRate',1, 'averageSpikeWaveform',12);


%% Once a set of neurons have been documented for each muscle, put this info
% together in a struct called 'common'



% for significance -- plot baseline samples per bin distribution, if normal
% then t-test, otherwise boostrap confidence interval. Also on PSTH plot,
% plot dashed line for baseline noise w/ conf interval (compute from beginning of file)

% in average spike waveform -- plot all other waveforms (mean, conf int) as
% well)

% randomely split file into two -- see if same neurons are in same file
