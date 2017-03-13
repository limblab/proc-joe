%script to set input data and execute data processing
%% process stimulation artifacts:
pwd = cd;
folderpath='D:\Lab\Data\StimArtifact\Chips\';
functionName='processStimArtifact';

inputData.task='tasknone';
inputData.ranBy='ranByTucker'; 
inputData.array1='arrayS1'; 
inputData.monkey='monkeyChips';
inputData.mapFile='mapFileD:\Lab\Data\MapFiles\Chips_Left_S1\SN 6251-001455.cmp'; % chips mapfile location

inputData.badChList=1:32;
inputData.interpulse=.000053;%in s
inputData.pWidth1=.0002;
inputData.pWidth2=.0002;

inputData.windowSize=30*20;%in points
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
outputData.artifactData = artifactData; clear artifactData;
outputData.chList = chList; clear chList;
outputData.eList = eList; clear eList;
outputData.posList = posList; clear posList;

%% add neuron like waveforms (triangles lol) to the artifact data
% add waves after the artifact -- data point inputData.presample and beyond
% waves should last roughly 0.5ms, sample rate is 30000
addWaveforms = 1;
sampleRate = 30000; % hz
timeWave = 0.5; % ms
minTimeWave = 0.2;
stdTimeWave = 0.1;
numWaves = 2;
outputData.artifactData.fakeWaves = cell(numel(outputData.eList),1);
outputData.artifactData.times = [];
if(addWaveforms)
    for i = 1:numel(outputData.eList) % for each electrode
        % build waveform
        lengthWave = ceil(max(minTimeWave,(timeWave+normrnd(0,stdTimeWave)))/1000*sampleRate);
        pointsDownDeflect = ceil(lengthWave*min(max(0.1,0.5+normrnd(0,0.1)),0.7));
        pointsUpDeflect = floor((lengthWave-pointsDownDeflect)/2);
        pointsToZero = max(1,lengthWave-pointsDownDeflect-pointsUpDeflect);
        fakeWave = zeros(pointsDownDeflect+pointsUpDeflect+pointsToZero-2,1);
        ampDown = -10;
        ampUp = 10;
        fakeWave(1:pointsDownDeflect) = linspace(0,ampDown,pointsDownDeflect);
        fakeWave(pointsDownDeflect:pointsDownDeflect+pointsUpDeflect-1) = linspace(ampDown,ampUp,pointsUpDeflect);
        fakeWave(pointsDownDeflect+pointsUpDeflect-1:pointsDownDeflect+pointsUpDeflect+pointsToZero-2) = linspace(ampUp,0,pointsToZero);
        for j = 1:size(outputData.artifactData.artifact,2) % for each stimulation
            times = floor((inputData.windowSize-4*length(fakeWave))*(rand(numWaves,1))+inputData.presample+length(fakeWave));
            for k = 1:numWaves
                fakeWaveRnd = fakeWave + normrnd(0,ampUp/25,size(fakeWave));
                outputData.artifactData.artifact(i,j,times(k):times(k)+length(fakeWave)-1) = squeeze(outputData.artifactData.artifact(i,j,times(k):times(k)+length(fakeWave)-1)) + fakeWaveRnd;
                outputData.artifactData.times(i,j,:) = times;
            end
        end
        outputData.artifactData.fakeWaves{i,:} = fakeWave';
    end
end
%% perform filtering step 
% filter combinations to use
highPassCutoff = [-1,500,2500,5000];
lowPassCutoff = [-1,1000,5000,10000];
filterOrder = [1];
% settling time matrix
settlingTimeMat = [];
filterParamsMat = [];
sampRate = 30000; % hz
for k = filterOrder
    for i = 1:numel(highPassCutoff)
        for j = 1:numel(lowPassCutoff)
            fhigh = highPassCutoff(i); % hz
            flow = lowPassCutoff(j); % hz
            doNothing = 0;
            noFilter = 0;
            filterStruct.userFilter = 1;
            if(flow==-1 && fhigh == -1) % do nothing because invalid filter
                doNothing = 0;
                noFilter = 1;
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

            if(~doNothing && ~noFilter) % filter and save png and .mat file
                outputDataFiltered = filterArtifactData(outputData,'filter',filterStruct);
                settlingTimeMat(end+1,:,:) = getSettlingTime(outputDataFiltered, inputData);
                % plot results and save png
                plotArtifactsAllStimChannels(outputDataFiltered,inputData,folderpath,'Name',filterStruct.name,'noPlots',1);
                save(strcat(cd,filesep,'Raw_Figures',filesep,filterStruct.name(2:end)),'outputDataFiltered','inputData');
                filterParamsMat(end+1,:) = [fhigh,flow,k];

            elseif(~doNothing && noFilter) % no filter, save .mat file
                settlingTimeMat(end+1,:,:) = getSettlingTime(outputData, inputData);
                plotArtifactsAllStimChannels(outputData,inputData,folderpath,'Name','_noFilter','noPlots',1);
                save(strcat(cd,filesep,'Raw_Figures',filesep,'noFilter'),'outputData','inputData');
                filterParamsMat(end+1,:) = [fhigh,flow,k];

            end
            
        end
    end
end

meanSettlingTime = mean(settlingTimeMat,3);

disp('Done plotting things');

%% uhh see if we can recover added waves?
% load file
load(strcat(cd,'\Raw_Figures\High_Order2_5000Hz.mat'));
load(strcat(cd,'\Raw_Figures\noFilter.mat'));
if(addWaveforms)
    figure();
    ch = 50;
    artNum = 10;
    plot(squeeze(outputData.artifactData.artifact(ch,artNum,:)),'k')
    hold on
    plot(squeeze(outputDataFiltered.artifactData.artifact(ch,artNum,:)),'r')
end