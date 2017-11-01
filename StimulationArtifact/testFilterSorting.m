%% process stimulation artifacts:
pwd = cd;
folderpath= 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20170621_3000stimuli\';
% folderpath='D:\Lab\Data\StimArtifact\Chips_one\';
functionName='processStimArtifact';

%%
warning('off')
inputData.task='taskRW';
inputData.ranBy='ranByTucker'; 
inputData.array1='arrayS1'; 
inputData.monkey='monkeyMihili';
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

%% run all for different positions of the waveform
% indexesToPlaceNeurons = [110,112,115,120,130];
indexesToPlaceNeurons = [80];
ampWave = 150/8;
tolerance = 3;
thresholdMult = 4;
preOffset = 33;
postOffset = 14; % these two should add up to 48-1?
removalSteps = {'','','','';...
    '0','High250_Order6','','';...
    '0','High500_Order6','',''};
% removalSteps = {'0','High500_Order6','',''};
removalSteps = buildRemovalAttemptName(removalSteps);

clear results;
clear outputDataFiltered;
for waveNum = 1:numel(indexesToPlaceNeurons)
    waveIdx = indexesToPlaceNeurons(waveNum);
    % load in artifact data generated above for filtering.
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
    
    % add neuron like waveforms to the artifact data
    % add waves after the artifact -- data point inputData.presample and beyond
    numWaves = 1;

%     outputData = addArtificialNeurons(outputData, neuronMeanWave, ampWave, waveIdx, numWaves);
    % perform artifact removal step
    outputDataFilteredTemp = removeArtifacts(outputData,inputData,removalSteps,folderpath);

    % threshold and store spike crossings
    nevData.ts = [];
    nevData.waveforms = [];
    nevData.elec = [];
    for tempIdx = 1:numel(outputDataFilteredTemp)
        for art = 1:numel(outputDataFilteredTemp{1,tempIdx}.artifactData)
            [tstemp,spikeWavestemp, chanstemp, lengthOut,shouldFind] = getThresholdCrossingData(outputDataFilteredTemp{1,tempIdx}.artifactData(art),...
                thresholdMult,preOffset,postOffset,waveIdx+14,art);
            % write nev file with above data
            nevData.ts = [nevData.ts; tstemp(1:lengthOut)];
            nevData.waveforms = [nevData.waveforms(:,:);spikeWavestemp(1:lengthOut,:)*8];
            nevData.elec = [nevData.elec;chanstemp(1:lengthOut)];
            packetWidth = (postOffset+preOffset+1)*2 + 8; % num bytes in data = length of data*2 + 8 bytes
            nevFilename = strcat('Han_20170516');
            mapFilename = inputData.mapFile(8:end);
            comments = '';
        end
        % sort spikes based on time
        [nevData.ts,sortTs] = sort(nevData.ts);
        nevData.waveforms = nevData.waveforms(sortTs,:);
        nevData.elec = nevData.elec(sortTs,:);
        disp('writing NEV file')
        tic
        writeNEV(nevData, packetWidth, nevFilename, mapFilename, comments)
        toc
    end

%     clear outputDataFilteredTemp;
end
disp('done writing NEV file')

%% make cds of sorted units from above
pwd=cd;
folderpath = 'D:\Lab\proc-joe\StimulationArtifact\';
filename = nevFilename;
labnum = 6;
monkey = 'monkeyHan';
ranBy = 'ranByTucker';
array = 'arrayLeftS1';n
task = 'taskRW';
cds=commonDataStructure();
cds.file2cds([folderpath,filename],ranBy,array,monkey,labnum,'ignoreJumps',task);
cd(pwd);

%% go through sorted units in cds and find waveform in raw data, plot
for spike = 1:1%size(cds.units,2)
%     if(cds.units(spike).ID~=255 && cds.units(spike).ID~=0) % if not unsorted or invalid
        for t = 1:1:5%size(cds.units(spike).spikes,1))
            ts = cds.units(spike).spikes{t,1};
            % parse ts for art: 1+(art-1)*10+(stim*0.01)+time
            art = floor(ts/10)+1;
            stim = floor((mod(ts,10)-1)*100)+1;
            idx = round((ts-1-(art-1)*10-((stim-1)*0.01))*30000);
            chan = cds.units(spike).chan;
            data = squeeze(outputData.artifactData(art).artifact(chan,stim,:));
            plotX = (0:1:(numel(data)-1))/30000*1000;
            figure;
            plot(plotX,data)
            hold on
            yline(1) = data(idx) + 1;
            yline(2) = data(idx) + 100;
            plot([plotX(idx),plotX(idx)],[yline(1),yline(2)],'r','linewidth',2);
        end
%     end 
end

%% writeNev testing
numSamples = 10000;
nevData.ts = linspace(0,60*60,numSamples);
nevData.elec = 1+zeros(numSamples,1);
nevData.waveforms = zeros(numSamples,48);
packetWidth = 48*2 + 8;
comments = '';
nevFilename = 'testingNev';
mapFilename = inputData.mapFile(8:end);
tic
writeNEV(nevData, packetWidth, nevFilename, mapFilename, comments)
toc
