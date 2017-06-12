function [ outputDataFiltered ] = recoverWaveformsPostFiltering(outputDataFiltered, inputData, thresholdMult, tolerance, filter, neuronMeanWave)

% initialize things
shift = computeShift(neuronMeanWave,filter, inputData);
electrodeStimCounter = 0;

artNum=numel(outputDataFiltered.artifactData);
neuronsFound(:,:,:) = zeros(artNum,size(outputDataFiltered.artifactData(1).artifact,1),size(outputDataFiltered.artifactData(1).artifact,2));
neuronsTotal(:,:,:) = zeros(artNum,size(outputDataFiltered.artifactData(1).artifact,1),size(outputDataFiltered.artifactData(1).artifact,2));
thresholdCrossings(:,:,:) = zeros(artNum,size(outputDataFiltered.artifactData(1).artifact,1),size(outputDataFiltered.artifactData(1).artifact,2));
neuronWaves(:,:,:,:) = zeros(artNum,size(outputDataFiltered.artifactData(1).artifact,1),size(outputDataFiltered.artifactData(1).artifact,2),...
    48); % waveforms have 48 datapoints, hard coding numbers is bad but whatever for right now
stimChanel = [];

% try to recover waves
tic
for art = 1:numel(outputDataFiltered.artifactData)
    stimChannel(art) = outputDataFiltered.artifactData(art).stimChannel;
    for i = 1:numel(outputDataFiltered.eList) 
        % check if bad channel
        if(sum(inputData.badChList-i == 0) == 0)
            for j = 1:size(outputDataFiltered.artifactData(art).artifact,2)
                % threshold based on rms (-4*rms)                       
                [total,numFound, thresholdCrossings] = thresholdData(squeeze(outputDataFiltered.artifactData(art).artifact(i,j,:)),...
                    outputDataFiltered.artifactData(art).fakeWaveTimes(i,j,:), abs(outputDataFiltered.artifactData(art).threshold(i,j))*thresholdMult,tolerance,shift,inputData);

                neuronsTotal(art,i,j) = total;
                neuronsFound(art,i,j) = numFound;
                fakeWaveTime = outputDataFiltered.artifactData(art).fakeWaveTimes(i,j)+shift; % shift by an amount based on the filter
                neuronWaves(art,i,j,:) = outputDataFiltered.artifactData(art).artifact(i,j,fakeWaveTime-14+shift:fakeWaveTime+33+shift);
                thresholdCrossings(art,i,j) = numel(thresholdCrossings);
                electrodeStimCounter = electrodeStimCounter+1;
            end
        end
    end
    percentFound(art) = sum(sum(neuronsFound(art,:,:)))/sum(sum(neuronsTotal(art,:,:)));
    meanThresholdCrossings(art) = sum(sum(thresholdCrossings(art,:,:)))/electrodeStimCounter;
end
toc
recoveredStruct.percentFound = percentFound;
recoveredStruct.meanThresholdCrossings = meanThresholdCrossings;
recoveredStruct.neuronsFound = neuronsFound;
recoveredStruct.neuronsTotal = neuronsTotal;
recoveredStruct.neuronWaves = neuronWaves;
recoveredStruct.stimChannel = stimChannel;


outputDataFiltered.recoveredWaveforms = recoveredStruct;

end

