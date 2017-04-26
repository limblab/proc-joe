%% process stimulation artifacts:
pwd = cd;
folderpath='D:\Lab\Data\StimArtifact\Mihili\lowGain_day2\noLineCancelling\chan69\';
% folderpath='D:\Lab\Data\StimArtifact\Chips_one\';
functionName='processStimArtifact';

%%
inputData.task='taskRW';
inputData.ranBy='ranByTucker'; 
inputData.array1='arrayPMD'; 
inputData.monkey='monkeyMihili';
inputData.mapFile='mapFileD:\Lab\Data\MapFiles\Mihili_Left_M1\Mihili Left M1 SN 1025-001452.cmp'; % chips mapfile location
% inputData.mapFile = 'mapFileD:\Lab\Data\MapFiles\Chips_Left_S1\SN 6251-001455.cmp';
inputData.badChList=0;
inputData.interpulse=.000053;%in s
inputData.pWidth1=.0002;
inputData.pWidth2=.0002;

inputData.windowSize=30*6;%in points
inputData.presample=50;%in points
inputData.plotRange=0.5;%in mV
inputData.lab=6;
inputData.useSyncLabel=[];

dataStruct2 = runDataProcessing(functionName,folderpath,inputData);
cd(pwd);


%% run all for different positions of the waveform
% indexesToPlaceNeurons = [110,112,115,120,130];
indexesToPlaceNeurons = [100];
ampWave = 0/8;
tolerance = 2;
thresholdMult = -2;
% removalSteps = {'','','','';...
%     '0','High250_Order6','','';...
%     '0','High500_Order6','',''};
removalSteps = {'','','',''};
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
    outputData = addArtificialNeurons(outputData, neuronMeanWave, ampWave, waveIdx, numWaves);
    
    % perform artifact removal step
    outputDataFilteredTemp = removeArtifacts(outputData,inputData,removalSteps,folderpath);
    
    % try to recover waves and store data
    for tempIdx = 1:numel(outputDataFilteredTemp)
        outputDataFiltered{tempIdx,waveNum} = ... 
            recoverWaveformsPostFiltering(outputDataFilteredTemp{1,tempIdx}, inputData, ...
            thresholdMult, tolerance, removalSteps{tempIdx,2},neuronMeanWave);
    end
    clear outputDataFilteredTemp;
end

%% plot percent found for each removal attempt
sampRate = 30000;
figure;
xData = [];
yData = [];
legendStr = {};
for removalAttempt = 1:size(outputDataFiltered,1)
    for waveIdx = 1:size(outputDataFiltered,2)
        xData(removalAttempt,waveIdx) = (indexesToPlaceNeurons(waveIdx)-inputData.presample+14)/sampRate;
        yData(removalAttempt,waveIdx) = mean(outputDataFiltered{removalAttempt,waveIdx}.recoveredWaveforms.percentFound);
    end
    plot(xData'*1000,yData','linewidth',2)
    legendStr{end+1,1} = removalSteps{removalAttempt,end};

end
ylim([0 1]);
xlabel('Time after stimulation (ms)');
ylabel('Percent found')
xlim([xData(1)*1000,xData(end)*1000])
l=legend(legendStr,'interpreter','none');
set(l,'interpreter','none')
set(l,'box','off')
formatForLee(gcf);

%% want to verify that neurons found are actually neurons....
% get all 48 data point windows around "found" neurons
eList = outputDataFiltered{1,1}.eList;
posList = outputDataFiltered{1,1}.posList;
chList = outputDataFiltered{1,1}.chList;
for waveIdx = 2:2%size(outputDataFiltered,2)
    for removalAttempt = 2:2%size(outputDataFiltered,1)
        for stimChan = 1:1%numel(outputDataFiltered{removalAttempt,waveIdx}.artifactData)
            % now we can go get all found neurons lol
            figure;
            for elec = 1:96
                posIdx=find(strcmp(eList,outputDataFiltered{removalAttempt,waveIdx}.artifactData(1).electrodeNames{elec}));
                eRow=posList(posIdx,1);
                eCol=posList(posIdx,2);
                subplot(10,10,(eRow-1)*10+eCol)
                hold on
                % plot all found waves for this electrode
                for stimulations = 1:1:size(outputDataFiltered{removalAttempt,waveIdx}.recoveredWaveforms.neuronsFound,3)
                    if(outputDataFiltered{removalAttempt,waveIdx}.recoveredWaveforms.neuronsFound(stimChan,elec,stimulations) > 0) % then plot wave
                        data = squeeze(outputDataFiltered{removalAttempt,waveIdx}.recoveredWaveforms.neuronWaves(stimChan,elec,stimulations,:));
                        data = data;
                        plot(data) 
                    end
                end
                ylim([-400/8 400/8])
             
            end
            suptitle('found')

            figure;
            for elec = 1:96
                posIdx=find(strcmp(eList,outputDataFiltered{removalAttempt,waveIdx}.artifactData(1).electrodeNames{elec}));
                eRow=posList(posIdx,1);
                eCol=posList(posIdx,2);
                subplot(10,10,(eRow-1)*10+eCol)
                hold on
                % plot all found waves for this electrode
                for stimulations = 1:1:size(outputDataFiltered{removalAttempt,waveIdx}.recoveredWaveforms.neuronsFound,3)
                    if(~(outputDataFiltered{removalAttempt,waveIdx}.recoveredWaveforms.neuronsFound(stimChan,elec,stimulations) > 0)) % then plot wave
                        data = squeeze(outputDataFiltered{removalAttempt,waveIdx}.recoveredWaveforms.neuronWaves(stimChan,elec,stimulations,:));
                        data = data - data(14);
                        plot(data) 
                    end
                end
                ylim([-400/8 400/8])
            end
            suptitle('missed')
        end
    end
end
% 
% %% plot percent found for each channel (results().foundStruct.neuronsFound)
% eList = outputDataFiltered.eList;
% posList = outputDataFiltered.posList;
% chList = outputDataFiltered.chList;
% for filter = 1:size(results(waveIdx).foundStruct.neuronWaves,1)
%     for stimChan = 1:size(results(1).foundStruct.neuronWaves,2)
%         % now we can go get all found neurons lol
%         % plot a 4x4 where each subplot is an electrode and all plots
%         % on that are the neurons found for that electrode. Do as many
%         % are necessary for this file.
%         figure;
%         for elec = 1:96
%             posIdx=find(strcmp(eList,outputDataFiltered.artifactData(1).electrodeNames{elec}));
%             eRow=posList(posIdx,1);
%             eCol=posList(posIdx,2);
%             h=subplot(10,10,10*(eRow-1)+eCol);
%             hold on
%             yData = zeros(numel(results),1);
%             xData = zeros(numel(results),1);
%             % plot percent found for each channel
%             for waveIdx = 1:numel(results)
%                 xData(waveIdx) = (results(waveIdx).waveIdx-100)/30000*1000;
%                 yData(waveIdx) = sum(squeeze(results(waveIdx).foundStruct.neuronsFound(filter,stimChan,elec,:)))/100;
%             end
%             f=plot(xData,yData,'linewidth',2,'color','k');
%             if(elec == results(waveIdx).foundStruct.stimChannel(stimChan))
%                 set(f,'color','r')
%             end
%             ylim([0 1])
%         end
%         suptitle(strrep(results(waveIdx).foundStruct.filenames{filter},'_',' '))
%     end
% end
% 
% %% compute and plot first two principle components for all found waveforms
% % as well as a set of noisy mean waveforms
% load('Neuron_data_pca.mat');
% load('unsortedWaves.mat');
% 
% foundNeuronData = zeros(50000,numel(neuronMeanWave));
% startWaveIdx = 1;
% filter = 1;
% numFound = 1;
% for waveIdx = startWaveIdx
%     for stimChan = 1:size(results(waveIdx).foundStruct.neuronWaves,2)
%         % now we can go get all found neurons lol
%         % plot a 4x4 where each subplot is an electrode and all plots
%         % on that are the neurons found for that electrode. Do as many
%         % are necessary for this file.
%         for elec = 1:96
%             % grab all found waves for this electrode
%             for stimulations = 1:size(results(waveIdx).foundStruct.neuronWaves,4)
%                 if(results(waveIdx).foundStruct.neuronWaves(filter,stimChan,elec,stimulations,1) ~= -1) % then plot wave
%                     data = squeeze(results(waveIdx).foundStruct.neuronWaves(filter,stimChan,elec,stimulations,:));
%                     data = data - mean(data);
%                     foundNeuronData(numFound,:) = data;
%                     numFound = numFound + 1;
%                     if(numFound == size(foundNeuronData,1))
%                         foundNeuronData = [foundNeuronData; zeros(1000,numel(neuronMeanWave))];
%                     end
%                 end   
%             end
%         end
%         
%     end
% end
% 
% numFoundDivisor = 2;
% foundNeuronData = foundNeuronData(1:floor(numFound/numFoundDivisor),:);
% numFound = floor(numFound/numFoundDivisor);
% numSimulated = 1000;
% simulatedNeuronData = repmat(neuronMeanWave,numSimulated,1) + randn(numSimulated,numel(neuronMeanWave)).*(2*repmat(neuronStdWave,numSimulated,1));
% unsortedWaves = unsortedWaves + randn(size(unsortedWaves))*mean(neuronStdWave);
% 
% markersize = 4;
% neuronData = [foundNeuronData; simulatedNeuronData; unsortedWaves];
% [coeff,score,latent]=pca(neuronData,'NumComponents',3);
% pcaNeuronData = neuronData*coeff;
% 
% plot(pcaNeuronData(numFound+numSimulated+1:end,1),pcaNeuronData(numFound+numSimulated+1:end,2),'.','markersize',markersize,'color','k')
% hold on
% plot(pcaNeuronData(numFound+1:numFound+numSimulated,1),pcaNeuronData(numFound+1:numFound+numSimulated,2),'.','markersize',markersize,'color','b')
% plot(pcaNeuronData(1:1:numFound,1),pcaNeuronData(1:1:numFound,2),'.','markersize',markersize,'color','r')
