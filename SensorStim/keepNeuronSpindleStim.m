function [keep] = keepNeuronSpindleStim(cds, neuronNumber, binEdges, binCounts, eventTimes, useRate, stimState, varargin)
% rejects or keeps neuron (nn in cds) based on the binEdges and counts
% provided

% keep neurons that have a "significant" difference in firing rate compared
% to all of the other bins. 

bootstrap = 1;
meanCounts = 1000;
plusMinus = [1000,1000];
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'bootstrap'
            bootstrap = varargin{i+1};
        case 'meanCounts'
            meanCounts = varargin{i+1};
        case 'plusMinus'
            plusMinus = varargin{i+1};
    end
end


binSize=binEdges(2)-binEdges(1);

if(bootstrap)
    [meanCounts, plusMinus] = bootstrapConfidenceInterval(cds,neuronNumber,eventTimes(1),binSize);
end

zeroIdx = find(binEdges==0);

keep = 0;
if(useRate)
    plusMinus=plusMinus/(binEdges(2)-binEdges(1));
end
threshold = 1.0;
stimSampRate = 1;
if(~isempty(cds.analog))
    stimSampRate = 1/(cds.analog{1,1}.t(2) - cds.analog{1,1}.t(1));
else
    stimSampRate = 1/(cds.lfp.t(2) - cds.lfp.t(1));
end
aveStimTime = (sum(stimState == 1))/numel(eventTimes)/stimSampRate;

if(sum(binCounts(zeroIdx+1:zeroIdx+floor(aveStimTime/binSize)) <= plusMinus(2)*threshold) <= 3 && ... % this is bootstrap constraint
        mean(binCounts(zeroIdx+1:zeroIdx+floor(aveStimTime/binSize))) > mean(binCounts(1:zeroIdx)) + 1.5*std(binCounts(1:zeroIdx)))
    keep = 1;
end

end

