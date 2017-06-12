% this script takes in post-filter waveforms and outputs potential pre-filter waveforms
% currently using fake data because of saving issues
load('Neuron_data_canonical')

[b,a] = butter(6,500/30000*2,'high');
neuronMeanWave = neuronMeanWave - mean(neuronMeanWave);
neuronWavesFiltered = fliplr(filter(b,a,fliplr(neuronWaves)')');

neuronMeanWaveFilt = mean(neuronWavesFiltered);
neuronStdWaveFilt = std(neuronWavesFiltered);

% we basically want to use a genetic algorithm to guess a bunch from random
% starting waveforms (48 point currently)

[minFitAll, countIterations,populationBest,populationAll] = runGA(neuronMeanWaveFilt,b,a);

%
figure;
subplot(2,1,1)
plot(neuronMeanWave,'k','linewidth',2);
hold on
plot(populationBest);
subplot(2,1,2)
plot(neuronMeanWaveFilt,'k','linewidth',1)
hold on
plot(fliplr(filter(b,a,fliplr(populationBest)')'))

%% see how many are close to the filtered waveform (index within 2 std)
lw = 2;

populationAllFilt = fliplr(filter(b,a,fliplr(squeeze(populationAll(:,:)))')');

score = sum(populationAllFilt-repmat(neuronMeanWaveFilt,size(populationAllFilt,1),1)<...
    repmat(2*neuronStdWaveFilt,size(populationAllFilt,1),1),2);

max(score)

wavesMask = (score == max(score));

wavesKeep = populationAll(wavesMask,:);
meanWavesKeep = mean(wavesKeep);
stdWavesKeep = std(wavesKeep);

plot(meanWavesKeep,'r','linewidth',lw)
hold on
plot(meanWavesKeep+stdWavesKeep*2,'r--','linewidth',lw)
plot(meanWavesKeep-stdWavesKeep*2,'r--','linewidth',lw)
plot(neuronMeanWave,'k','linewidth',lw)



