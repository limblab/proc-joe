%script to set input data and execute data processing
%% process stimulation artifacts:
pwd = cd;
folderpath='D:\Lab\Data\StimArtifact\Chips\';
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

inputData.windowSize=30*10;%in points
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
inputData = temp;
%% add neuron like waveforms (triangles lol) to the artifact data
% add waves after the artifact -- data point inputData.presample and beyond
% waves should last roughly 0.5ms, sample rate is 30000
load('Neuron_data.mat');
lambda = 0.6; % lambda for poisson distribution
addWaveforms = 1;
sampleRate = 30000; % hz

if(addWaveforms)
    for a = 1:numel(outputData.artifactData)
        outputData.artifactData(a).fakeWaveTimes = zeros(numel(outputData.eList),size(outputData.artifactData(a).artifact,20))-1;
        for i = 1:numel(outputData.eList) % for each electrode        
            for j = 1:size(outputData.artifactData(a).artifact,2) % for each stimulation
                % determine number of waves
                numWaves = poissrnd(lambda);
                times = floor((inputData.windowSize-inputData.presample-1-length(neuronMeanWave))*(rand(numWaves,1))+inputData.presample);
                for k = 1:numWaves % for each wave to add
                    wave = neuronMeanWave/3;
                    for w = 1:numel(wave) % for each point in wave, add noise
                        wave(w) = wave(w) + randn()*neuronStdWave(w)/3;
                    end
                    outputData.artifactData(a).artifact(i,j,times(k):times(k)+length(wave)-1) = squeeze(outputData.artifactData(a).artifact(i,j,times(k):times(k)+length(wave)-1))'+wave;
                    outputData.artifactData(a).fakeWaveTimes(i,j,1:numel(times)) = times;
                end
            end
        end
    end
end
%% perform filtering step 
% filter combinations to use
highPassCutoff = [1500];
lowPassCutoff = [2000,2500,3000];
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
                save(strcat(folderpath,filesep,'Raw_Figures',filesep,filterStruct.name(2:end)),'outputDataFiltered','inputData');
                filterParamsMat(end+1,:) = [fhigh,flow,k];

            elseif(~doNothing && noFilter) % no filter, save .mat file
                settlingTimeMat(end+1,:,:) = getSettlingTime(outputData, inputData);
                plotArtifactsAllStimChannels(outputData,inputData,folderpath,'Name','_noFilter','noPlots',1);
                save(strcat(folderpath,filesep,'Raw_Figures',filesep,'noFilter'),'outputData','inputData');
                filterParamsMat(end+1,:) = [fhigh,flow,k];

            end
            
        end
    end
end

meanSettlingTime = mean(settlingTimeMat,3);
meanSettlingTime = mean(meanSettlingTime(:,33:end),2);
disp('Done plotting things');



%% See if waves can be recovered?
thresholdMult = -4;
offset = 2;


filteredDataFilepath = strcat(folderpath,'Raw_Figures',filesep);
% get all file names
namesStruct = dir(strcat(filteredDataFilepath,'*.mat'));
filteredDataFilenames = cell(1,numel(namesStruct));
for i = 1:numel(namesStruct)
    filteredDataFilenames{i} = namesStruct(i).name;
end

foundStruct.neuronsFound = [];
foundStruct.neuronsTotal = [];
foundStruct.percentFound = [];
foundStruct.meanThresholdCrossings = [];
foundStruct.thresholdCrossings = [];
foundStruct.filenames = {};
foundStruct.filepath = filteredDataFilepath;
foundStruct.neuronMeanWave = neuronMeanWave;
foundStruct.neuronStdWave = neuronStdWave;
% load in each file, then see if waves were added
for fname = filteredDataFilenames % fname is the filename then
    load(strcat(filteredDataFilepath,fname{1}));
    if(isfield(outputDataFiltered.artifactData,'fakeWaveTimes'))
        for art = 1:numel(outputDataFiltered.artifactData)
            foundStruct.filenames{end+1} = fname{1};
            foundStruct.neuronsFound(end+1,art,:,:) = zeros(size(outputDataFiltered.artifactData(art).artifact,1),size(outputDataFiltered.artifactData(art).artifact,2));
            foundStruct.neuronsTotal(end+1,art,:,:) = zeros(size(outputDataFiltered.artifactData(art).artifact,1),size(outputDataFiltered.artifactData(art).artifact,2));
            foundStruct.thresholdCrossings(end+1,art,:,:) = zeros(size(outputDataFiltered.artifactData(art).artifact,1),size(outputDataFiltered.artifactData(art).artifact,2));
            counter = 0;

            outputDataFiltered.artifactData(art).foundWave = outputDataFiltered.artifactData(art).fakeWaveTimes*0; % initialize found neuron matrix
            for i = 33:numel(outputDataFiltered.eList) % 1:32 are bad yo
                for j = 1:size(outputDataFiltered.artifactData(art).artifact,2)
                    % threshold based on rms (-4*rms is what Tucker usually uses???)
                    artifact = outputDataFiltered.artifactData(art).artifact;
                    threshold = thresholdMult*rms(artifact(i,j,inputData.presample + inputData.windowSize/2:end));
                    for t = 1:size(outputDataFiltered.artifactData(art).fakeWaveTimes,3)
                        neuronIdx = outputDataFiltered.artifactData(art).fakeWaveTimes(i,j,t);
                        if(neuronIdx <= 0) % no neuron
                            outputDataFiltered.artifactData(art).foundWave(i,j,t) = -1;
                        elseif(sum(artifact(i,j,neuronIdx-offset:neuronIdx+offset) < threshold) > 0) % neuron crossed threshold
                            outputDataFiltered.artifactData(art).foundWave(i,j,t) = 1;    
                            foundStruct.neuronsTotal(end,art,i,j) = foundStruct.neuronsTotal(end,art,i,j) + 1;
                            foundStruct.neuronsFound(end,art,i,j) = foundStruct.neuronsFound(end,art,i,j) + 1;
                        else
                            foundStruct.neuronsTotal(end,art,i,j) = foundStruct.neuronsTotal(end,art,i,j) + 1;
                        end
                        foundStruct.thresholdCrossings(end,art,i,j) = sum(artifact(i,j,inputData.presample:end) < threshold);
                    end
                    counter = counter+1;
                end
            end
            if(art==1)
                foundStruct.percentFound(art,end+1) = sum(sum(foundStruct.neuronsFound(end,art,:,:)))/sum(sum(foundStruct.neuronsTotal(end,art,:,:)));
                foundStruct.meanThresholdCrossings(art,end+1) = sum(sum(foundStruct.thresholdCrossings(end,art,:,:)))/counter;
            else
                foundStruct.percentFound(art,end) = sum(sum(foundStruct.neuronsFound(end,art,:,:)))/sum(sum(foundStruct.neuronsTotal(end,art,:,:)));
                foundStruct.meanThresholdCrossings(art,end) = sum(sum(foundStruct.thresholdCrossings(end,art,:,:)))/counter;
            end
            save(strcat(filteredDataFilepath,fname{1}),'outputDataFiltered','inputData');
        end
    end
end
    
