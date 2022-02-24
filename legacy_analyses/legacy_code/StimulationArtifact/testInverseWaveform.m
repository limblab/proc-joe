% this script takes in post-filter waveforms and outputs potential pre-filter waveforms
% currently using fake data because of saving issues
load('Neuron_data_canonical')

[b,a] = butter(6,500/30000*2,'high');
neuronWavesFiltered = fliplr(filter(b,a,fliplr(neuronWaves)')');

neuronMeanWave = mean(neuronWavesFiltered);
neuronStdWave = std(neuronWavesFiltered);

% we basically want to guess a ton of possible waveforms, filter and check
% bounds -- then move those guesses towards the right thing -- Genetic Alg
numPop = 1000;
numRuns = 100;
pop = zeros(numPop,numRuns+1,48);
% init pop
rangeGuess = max(neuronMeanWave) - min(neuronMeanWave) + 10;
meanGuess = mean(neuronMeanWave);
pop(:,1,:) = rand(numPop,48)*rangeGuess - rangeGuess/2 + meanGuess;

for r = 1:numRuns
    % eval pop -- distance to neuronMeanWave 
    popFilt = fliplr(filter(b,a,fliplr(squeeze(pop(:,r,:)))')');
    scores = sum((popFilt-repmat(neuronMeanWave,numPop,1)).^2,2);
    % resample pop based on scores
    prob = scores/sum(scores);
    for mem = 1:numPop
        mem1 = ceil(rand()*numPop);
        mem2 = ceil(rand()*numPop);
        if(prob(mem1) > prob(mem2))
            pop(mem,r+1,:) = squeeze(pop(mem1,r,:));
        else
            pop(mem,r+1,:) = squeeze(pop(mem2,r,:));
        end
    end
    
    % crossover, random mutate
    for cross = 1:numPop*0.8
        mem1 = ceil(rand()*numPop);
        mem2 = ceil(rand()*numPop);
        
        pop(mem1,r+1,1:24) = squeeze(pop(mem2,r+1,1:24)) + rand(24,1)*200 - 100;
        pop(mem1,r+1,25:end) = squeeze(pop(mem2,r+1,25:end)) + rand(24,1)*200 - 100;
    end
end