function [ minFitAll, countIterations,populationWavesBest,populationWavesAll,N ] = runGA(neuronMeanWaveFiltered, b, a)

%% fit the data in trainingDataTrimmed
% A,B,C,D
% restrict a,b,c,d [-10,10]
% f(x) = 
global N maxBit popSize numPoints minRange maxRange data numElites numRandom bFilter aFilter

bFilter = b;
aFilter = a;
numElites = 0;
numRandom = ceil(popSize*0.1);
% generate data
data = neuronMeanWaveFiltered;
% each chromosome will have 6 N length bitstrings
N = 48; % will map the bit range to [min(data), max(data)];
numPoints = numel(neuronMeanWaveFiltered);
maxBit = 2^N - 1;
popSize = 250;
minRange = min(data)*2;
maxRange = max(data)*2;
threshold = 0;
maxIterations = 750;

minFitAll = zeros(maxIterations,1)-1;
%% Initialize Population
population = num2str(round(rand(popSize,N*numPoints)),'%d');
populationAll(:,:,1) = population;
%% Evaluate Population
allFit = evaluate(population); % returns fitness -> f(x) value, so lower fitness is better
[minFit,bestIdx] = min(allFit);
populationBest = population(bestIdx,:);
%% While loop
countIterations = 0;
while(minFit > threshold && countIterations < maxIterations)
    intPop = selectIntermediatePopulation(population, allFit);
    intPop = crossoverPopulation(intPop);
    population = mutatePopulation(intPop);
    allFit = evaluate(population);
    [minFit,bestIdx] = min(allFit);
    countIterations = countIterations+1;
    minFitAll(countIterations) = minFit;
    populationBest = population(bestIdx,:);
    populationAll(:,:,countIterations) = population;

end

populationWavesBest = zeros(1,numPoints);
populationWavesAll = zeros(popSize,numPoints,maxIterations);
for i = 1:numPoints
    idxStart = 1+(i-1)*N;
    idxEnd = N*i;
    populationWavesBest(1,i) = bin2dec(populationBest(:,idxStart:idxEnd))/maxBit*(maxRange-minRange) + minRange;
    for s = 1:maxIterations
        populationWavesAll(:,i,s) = bin2dec(populationAll(:,idxStart:idxEnd,s))/maxBit*(maxRange-minRange) + minRange;
    end
end

populationWavesBest = populationWavesBest - mean(populationWavesBest);

ptemp = zeros(size(populationWavesAll,1)*size(populationWavesAll,3),size(populationWavesAll,2));
idxStart = 1;
for iter = 1:size(populationWavesAll,3)
    idxEnd = idxStart+size(populationWavesAll,1)-1;
    ptemp(idxStart:idxEnd,:) = squeeze(populationWavesAll(:,:,iter));
    idxStart = idxEnd + 1;
end

populationWavesAll = ptemp;
populationWavesAll = populationWavesAll - mean(populationWavesAll,2);
end

