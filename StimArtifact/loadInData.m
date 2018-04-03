folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\testingCode\';

inputData.monkey = 'monkeyHan';
inputData.array =  'arrayLeftS1';
inputData.labnum = 6;
inputData.task = 'taskCObump';
inputData.ranBy = 'ranByJoseph';
inputData.mapFile = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

cd(folderpath)
fileList = dir('*spikesExtracted.nev*');
outputDataFileList = dir('*outputData*');


%%
fileNum = 2;
disp(fileList(fileNum).name)
cds = commonDataStructure();
cds.file2cds([folderpath fileList(fileNum).name],inputData.monkey,inputData.array,...
    inputData.labnum,inputData.task,inputData.ranBy,inputData.mapFile);

load(outputDataFileList(fileNum).name);

cdsTemp = [];
cdsTemp.meta = cds.meta;
cdsTemp.kin = cds.kin;
cdsTemp.force = cds.force;
cdsTemp.lfp = cds.lfp;
cdsTemp.emg = cds.emg;
cdsTemp.analog = cds.analog;
cdsTemp.triggers = cds.triggers;
cdsTemp.units = cds.units;
cdsTemp.trials = cds.trials;
cdsTemp.aliasList = cds.aliasList;
cdsTemp.operationLog = cds.operationLog;
stimOnIdx = find(diff(cds.analog{1}.ainp16-mean(cds.analog{1}.ainp16)>3)>.5);
cdsTemp.stimOn = cds.analog{1}.t(stimOnIdx);
% cdsTemp.stimOn = (outputData.stimInfo.stimOn/30000)';
% cdsTemp.stimOff = (outputData.stimInfo.stimOff/30000)';
cdsTemp.waveforms = outputData.waveforms;
cdsTemp.artifactData = outputData.artifactData;
cds = cdsTemp;