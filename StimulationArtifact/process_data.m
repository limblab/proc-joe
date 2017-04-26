%% process stimulation artifacts:
pwd = cd;
folderpath='D:\Lab\Data\StimArtifact\Chips_one\';
functionName='processStimArtifact';

%%
inputData.task='tasknone';
inputData.ranBy='ranByTucker'; 
inputData.array1='arrayS1'; 
inputData.monkey='monkeyChips';
inputData.mapFile='mapFileD:\Lab\Data\MapFiles\Chips_Left_S1\SN 6251-001455.cmp'; % chips mapfile location

inputData.badChList=1:32;
inputData.interpulse=.000053;%in s
inputData.pWidth1=.0002;
inputData.pWidth2=.0002;

inputData.windowSize=30*6;%in points
inputData.presample=100;%in points
inputData.plotRange=0.2;%in mV
inputData.lab=6;
inputData.useSyncLabel=[];

dataStruct2 = runDataProcessing(functionName,folderpath,inputData);
cd(pwd);

%% load in artifact data generated above for filtering.
load(strcat(folderpath,'Input_Data\','Input_structure.mat'));
load(strcat(folderpath,'Output_Data\','artifactData.mat'));
load(strcat(folderpath,'Output_Data\','chList.mat'));
load(strcat(folderpath,'Output_Data\','eList.mat'));
load(strcat(folderpath,'Output_Data\','posList.mat'));
load('Neuron_data_canonical.mat');
outputData.artifactData = artifactData; clear artifactData;
outputData.chList = chList; clear chList;
outputData.eList = eList; clear eList;
outputData.posList = posList; clear posList;
inputData = temp;
%% add neuron like waveforms to the artifact data
% add waves after the artifact -- data point inputData.presample and beyond
% waves should last roughly 0.5ms, sample rate is 30000
ampWave = 200;
waveIdx = 200;
numWaves = 1;
outputData = addArtificialNeurons(outputData, neuronMeanWave, ampWave, waveIdx, numWaves);

%% Template Subtraction
outputData = performTemplateSubtraction(outputData);

%% perform filtering step 
disp('Start Filtering Step');
% filter combinations to use
highPassCutoff = [-1,200,400,600];
lowPassCutoff = [-1,1300,1500,1700,2000];
filterOrder = [6,8];
[outputDataFiltered, filterResults] = performFilteringStep(outputData, inputData, highPassCutoff, lowPassCutoff, filterOrder, folderpath);
disp('End filtering step');
disp('')


%% See if waves can be recovered?
disp('Start wave recovery');

tolerance = 2;
thresholdMult = -4;
filteredDataFilepath = strcat(folderpath,'Raw_Figures',filesep);
neuronsFound = recoverWaveformsPostFiltering(outputDataFiltered, inputData, filteredDataFilepath, thresholdMult, neuronMeanWave, tolerance);

disp(num2str(rand()*20))
disp('End wave recovery');

