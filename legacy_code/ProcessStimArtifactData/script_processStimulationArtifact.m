%% process stimulation artifacts:
pwd = cd;
% folderpath='C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20180310_chic201802\';
folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20180311_stimswitch\';
functionName='processStimArtifact_filter'; % includes acausal filtered data

%%
inputData.saveFigures = 0;
inputData.dukeBoardChannel = -1;
inputData.dukeBoardLabel = 'ainp15';

inputData.task='taskCObump';
inputData.ranBy='ranByTucker'; 
inputData.array1='arrayS1'; 
inputData.monkey='monkeyHan';
inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp'; % mapfile location

inputData.badChList=[];
inputData.interpulse=.000053;%in s
inputData.pWidth1=.0002;
inputData.pWidth2=.0002;

inputData.windowSize=30*10;%in points
inputData.presample=100;%in points
inputData.plotRange=8;%in mV
inputData.plotRangeFiltered=0.5;%in mV
inputData.lab=6;
inputData.useSyncLabel=[];

[output_data] = runDataProcessing(functionName,folderpath,inputData);
disp('done processing')
cd(pwd);


