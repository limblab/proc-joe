% this script takes in post-filter waveforms and outputs potential pre-filter waveforms
% currently using fake data because of saving issues
load('Neuron_data_canonical')

[b,a] = butter(6,500/30000*2,'high');
neuronMeanWave = neuronMeanWave - mean(neuronMeanWave);
neuronWavesFiltered = fliplr(filter(b,a,fliplr(neuronWaves)')');

neuronMeanWaveFilt = mean(neuronWavesFiltered);
neuronStdWaveFilt = std(neuronWavesFiltered);

% we basically want to use grad descent
alpha = 0.01;
beta = 0.1;
% [waves, numWaves] = runGradDescent(neuronMeanWaveFilt,b,a,'lineSearch',[alpha,beta]);
tic
[waves, numWaves] = runGradDescent(neuronMeanWaveFilt,b,a);
toc
%
figure;
subplot(2,1,1)
plot(neuronMeanWave,'k','linewidth',2);
hold on
plot(waves(numWaves,:),'r');
subplot(2,1,2)
plot(neuronMeanWaveFilt,'k','linewidth',1)
hold on
plot(fliplr(filter(b,a,fliplr(waves(numWaves,:))')'),'r')

%% see how many are close to the filtered waveform (index within 2 std)
lw = 2;

wavesFilt = fliplr(filter(b,a,fliplr(squeeze(waves(:,:)))')');

score = sum(wavesFilt-repmat(neuronMeanWaveFilt,size(wavesFilt,1),1)<...
    repmat(1*neuronStdWaveFilt,size(wavesFilt,1),1),2);

max(score)

wavesMask = (score >= max(score));

wavesKeep = waves(wavesMask,:);
meanWavesKeep = mean(wavesKeep);
stdWavesKeep = std(wavesKeep);
figure
plot(meanWavesKeep,'r','linewidth',lw)
hold on
plot(meanWavesKeep+stdWavesKeep*2,'r--','linewidth',lw)
plot(meanWavesKeep-stdWavesKeep*2,'r--','linewidth',lw)
plot(neuronMeanWave,'k','linewidth',lw)



