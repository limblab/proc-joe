%% process stimulation artifacts:
pwd = cd;
folderpath='C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20180311_stimswitch\';
functionName='processStimArtifact';

%%
inputData.task='taskCObump';
inputData.ranBy='ranByTucker'; 
inputData.array1='arrayS1'; 
inputData.monkey='monkeyHan';
inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp'; % chips mapfile location

inputData.badChList=[];
inputData.interpulse=.000053;%in s
inputData.pWidth1=.0002;
inputData.pWidth2=.0002;

inputData.windowSize=30*8;%in points
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
highPassCutoff = [-1];
lowPassCutoff = [2000];
filterOrder = [6];
[outputDataFiltered] = performFilteringStep(outputData, inputData, highPassCutoff, lowPassCutoff, filterOrder, folderpath);
disp('End filtering step');
disp('')

%% plot artifact data before and after filtering
% for chan = 1:96
chan = 12;
    fileprefix = 'StimTesting_20170625_chan3stim_chan';
    filename = strcat(fileprefix,num2str(chan),'record');
    stimuli = 1:1:10;
    xData = ((1:size(outputData.artifactData.artifact,3))-101)/30;
    yData = squeeze(outputData.artifactData.artifact(chan,stimuli,:));
    f=figure();
    plot(xData,yData)
    title('No Filter')
    xlabel('Time after stimulation onset (ms)');
    ylabel('Voltage (\muV)')
    xlim([-1,3]);
    formatForLee(gcf);
%     saveFiguresLIB(gcf,'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\stimTesting\20170625\',filename);
%     ylim([-500,500])
%     saveFiguresLIB(gcf,'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\stimTesting\20170625\',strcat(filename,'_zoomed'));

    filtername = '6thOrder_2000HzHighPass';
    [b,a]=butter(6,2000/(30000/2),'low');
    yDataFiltered = filtfilt(b,a,yData')';
    figure();
    plot(xData,yDataFiltered)
    title('6th Order, 2kHz low pass filtered')
    xlabel('Time after stimulation onset (ms)');
    ylabel('Voltage (\muV)')
    xlim([-1,3]);
    formatForLee(gcf);
    saveFiguresLIB(gcf,'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\stimTesting\20170625\',strcat(filename,'_',filtername));
%     ylim([-500,500])
%     saveFiguresLIB(gcf,'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\stimTesting\20170625\',strcat(filename,'_',filtername,'_zoomed'));
%     close all
% end
%% See if waves can be recovered?
disp('Start wave recovery');

tolerance = 2;
thresholdMult = -4;
filteredDataFilepath = strcat(folderpath,'Raw_Figures',filesep);
neuronsFound = recoverWaveformsPostFiltering(outputDataFiltered, inputData, filteredDataFilepath, thresholdMult, neuronMeanWave, tolerance);

disp(num2str(rand()*20))
disp('End wave recovery');

%% plot neuron waves and filtered version
filename = 'Neuron_';
nn = 99;

waves = cds.units(nn).spikes{:,2:end};
plot(waves(1:10,:)')
rawIdx = zeros(size(waves,1),1);
for r = 1:numel(rawIdx)
    rawIdx(r) = getRawDataIdx(cds.units(nn).spikes{r,1},cds.units(nn).chan,cds.rawData.ts,cds.rawData.elec);
end

%%
rawIdx(rawIdx==-1) = [];
%%
    xData = ((1:size(cds.rawData.waveforms,2))-1)/30;
    yData = squeeze(cds.rawData.waveforms(rawIdx(1:50),:));
    f=figure();
    plot(xData,(yData'-mean(yData')))
    ylim([-200,200])
    %%
    
    title('No Filter')
    xlabel('Time after stimulation onset (ms)');
    ylabel('Voltage (\muV)')
    xlim([-1,3]);
    formatForLee(gcf);
    saveFiguresLIB(gcf,'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\stimTesting\20170625\',filename);
    ylim([-500,500])
    saveFiguresLIB(gcf,'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\stimTesting\20170625\',strcat(filename,'_zoomed'));

    filtername = '6thOrder_2000HzHighPass';
    [b,a]=butter(6,2000/(30000/2),'high');
    yDataFiltered = filtfilt(b,a,yData')';
    figure();
    plot(xData,yDataFiltered)
    title('6th Order, 2kHz high pass filtered')
    xlabel('Time after stimulation onset (ms)');
    ylabel('Voltage (\muV)')
    xlim([-1,3]);
    formatForLee(gcf);
    saveFiguresLIB(gcf,'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\stimTesting\20170625\',strcat(filename,'_',filtername));
    ylim([-500,500])
    saveFiguresLIB(gcf,'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\stimTesting\20170625\',strcat(filename,'_',filtername,'_zoomed'));
    close all