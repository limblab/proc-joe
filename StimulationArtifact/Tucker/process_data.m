%script to set input data and execute data processing
%% process stimulation artifacts:
folderpath='D:\Lab\Data\StimArtifact\Chips\';
functionName='processStimArtifact';

inputData.task='tasknone';
inputData.ranBy='ranByTucker'; 
inputData.array1='arrayS1'; 
inputData.monkey='monkeyChips';
inputData.mapFile='mapFileD:\Lab\Data\MapFiles\Chips_Left_S1\SN 6251-001455.cmp'; % chips mapfile location

inputData.windowSize=30*10;%in points
inputData.presample=5;%in points
inputData.plotRange=0.5;%in mV
inputData.lab=6;
inputData.useSyncLabel=[];

dataStruct2 = runDataProcessing(functionName,folderpath,inputData);

%% load in artifact data generated above for filtering.
load(strcat(folderpath,'Input_Data\','Input_structure.mat'));
load(strcat(folderpath,'Output_Data\','artifactData.mat'));
load(strcat(folderpath,'Output_Data\','chList.mat'));
load(strcat(folderpath,'Output_Data\','eList.mat'));
load(strcat(folderpath,'Output_Data\','posList.mat'));
outputData.artifactData = artifactData; clear artifactData;
outputData.chList = chList; clear chList;
outputData.eList = eList; clear eList;
outputData.posList = posList; clear posList;

% perform filtering step
highPassCutoff = [-1,5,];%10,50,100,250,1000];
lowPassCutoff = [-1,100];%,1000,2000,5000];
filterOrder = [1,2];%,3];
sampRate = 30000; % hz
for k = filterOrder
    for i = 1:numel(highPassCutoff)
        for j = 1:numel(lowPassCutoff)
            fhigh = highPassCutoff(i); % hz
            flow = highPassCutoff(j); % hz
            doNothing = 0;
            filterStruct.userFilter = 1;
            if(flow==-1 && fhigh == -1) % do nothing because invalid filter
                doNothing = 1;
            elseif(fhigh == -1) % only use low pass
                [filterStruct.b,filterStruct.a] = butter(2,flow/(sampRate/2),'low');
                filterStruct.name = strcat('_Low_Order',num2str(k),'_',num2str(flow),'Hz');
            elseif(flow == -1) % only use high pass
                [filterStruct.b,filterStruct.a] = butter(2,fhigh/(sampRate/2),'high');
                filterStruct.name = strcat('_High_Order',num2str(k),'_',num2str(fhigh),'Hz');
            elseif(flow <= fhigh) % do nothing because invalid filter
                doNothing = 1;
            else % band pass
                [filterStruct.b,filterStruct.a] = butter(2,[fhigh,flow]/(sampRate/2),'bandpass');
                filterStruct.name = strcat('_Bandpass_Order',num2str(k),'_',num2str(fhigh),'Hz',num2str(flow),'Hz');
            end

            if(~doNothing) % filter and save png
                outputDataFiltered = filterArtifactData(outputData,'filter',filterStruct);
                % plot results and save png
                plotArtifactsAllStimChannels(outputDataFiltered,inputData,folderpath,'Name',filterStruct.name);
                close all;
            end
        end
    end
end

