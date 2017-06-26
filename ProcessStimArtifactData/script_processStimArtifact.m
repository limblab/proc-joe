%% process stimulation artifacts:
pwd = cd;
folderpath= 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_processed\20170616\';
functionName='processStimArtifactData';

%% generates _cds and _crossings files
warning('off')
inputData.task='taskCO';
inputData.ranBy='ranByJoseph'; 
inputData.array1='arrayS1'; 
inputData.monkey='monkeyHan';
inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp'; % chips mapfile location
inputData.badChList=0;
inputData.interpulse=.000053;%in s
inputData.pWidth1=.0002;
inputData.pWidth2=.0002;

inputData.windowSize=30*10;%in points
inputData.presample=1;%in points
inputData.plotRange=0.5;%in mV
inputData.lab=6;
inputData.useSyncLabel=[];

dataStruct2 = runDataProcessing(functionName,folderpath,inputData);
cd(pwd);
warning('on')

%% blank artifact and writeNev from _nevData file
pwd = cd;
cd(folderpath);
fileList=dir('*_nevData*');
load(fileList(1).name);
fileList=dir('*_cds*');
load(fileList(1).name);
packetWidth = 104;
filename = strcat(fileList(1).name(1:end-14),'_merged_testing');
mapFilename = inputData.mapFile(8:end);
comments = '';
[nevData] = blankStimulationArtifact(cds,nevData);
writeNEV(nevData, packetWidth, filename, mapFilename, comments )
cd(pwd);
%% sort *_merged and call it *_merged-s

%% load in *_merged-s, save the units and move that to cds 
pwd = cd;
cd(folderpath);
fileList=dir('*cds*');
load(fileList(1).name);
fileList = dir('*_merged-s.NEV*');
filename = fileList(1).name;
labnum = 6;
monkey = 'monkeyHan';
ranBy = 'ranByJoseph';
array = 'arrayLeftS1';
task = 'taskRW';
spikes = commonDataStructure();
spikes.file2cds([folderpath,filename],ranBy,array,monkey,labnum,'ignoreJumps',task);
cds.units = spikes.units;
clear spikes
cd(folderpath)
save(strcat(filename(1:end-10),'_processed'),'cds','-v7.3');
cd(pwd);

%% copy relevant file by hand to folder on local computer for further analysis