function [keep] = keepNeuronGTOstim(cds, neuronNumber, binEdges, binCounts, stimStart, useRate, varargin)
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
    [meanCounts, plusMinus] = bootstrapConfidenceInterval(cds,neuronNumber,stimStart,binSize);
end

zeroIdx = find(binEdges==0);

keep = 0;
if(useRate)
    plusMinus=plusMinus/(binEdges(2)-binEdges(1));
end
threshold = 1.0;

% keep if count in bin is large enough and this is the max bin in the trial
% time
if(binCounts(zeroIdx+1) > plusMinus(2)*threshold && binCounts(zeroIdx+1)>mean(binCounts)+2*std(binCounts))
    keep = 1;
% elseif(binCounts(zeroIdx+2) > plusMinus(2)*threshold && max(binCounts) == binCounts(zeroIdx+2))
%     keep = 1;
end

end

