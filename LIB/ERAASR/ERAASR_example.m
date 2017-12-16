%% this script provides an example of ERAASR for our artifact data

% load in dataset
folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Chips_20171201';
mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';

pwd=cd;
cd(folderpath)
fileList = dir('*all_processed*');
fileNumber = 1;
load(fileList(fileNumber).name);
cd(pwd)

%% format artifact data into C x Time x Pulses x R and remove excess data
artMult = numel(cds.stimOn)/size(cds.artifactData.artifact,1);
artifactsKeep = cds.waveforms.waveSent(1:artMult:end) == 1;
artifact = cds.artifactData.artifact;
artifact = artifact(artifactsKeep,:,:);
% clear cds

artifact = permute(artifact,[2,3,4,1]);

%% setup opts for ERAASR

opts.K = [25,-1,10]; % number of PCs to keep for each step
opts.LAMBDA = [0,0,0]; % number of adjacent channels to remove
opts.BETA = [0,1]; % true or false
opts.ALIGN_TRIALS = 1;
%% run ERAASR

artifactCleaned = ERAASR(artifact,opts);

%% plot artifact before and after cleaning for a couple of channels
chanToPlot = 52;
figure
plot(squeeze(artifactCleaned.input(chanToPlot,:,1,1:10:100)))
figure
plot(squeeze(artifactCleaned.postTrials(chanToPlot,:,1,1:10:100)))













