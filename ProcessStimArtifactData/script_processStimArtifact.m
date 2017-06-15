%% process stimulation artifacts:
pwd = cd;
folderpath= 'R:\data\Han_13B1\Raw\bumpstim\20170614\';
% folderpath='D:\Lab\Data\StimArtifact\Chips_one\';
functionName='processStimArtifactData';

%%
warning('off')
inputData.task='taskCO';
inputData.ranBy='ranByJoseph'; 
inputData.array1='arrayS1'; 
inputData.monkey='monkeyHan';
inputData.mapFile='mapFileD:\Lab\Data\MapFiles\Han_Left_S1\SN 6251-001459.cmp'; % chips mapfile location
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