function [keep] = keepNeuronGTOstim(cds, nn, binEdges, binCounts)
% rejects or keeps neuron (nn in cds) based on the binEdges and counts
% provided

% keep neurons that have a "significant" difference in firing rate compared
% to all of the other bins. 

keep=1;
threshold = 2;

meanRate = mean(binCounts);
stdRate = std(binCounts);

idxBin = find(binEdges==0)+1;

if(binCounts(idxBin)>meanRate+stdRate*threshold)
    keep=1;
else
    keep=0;
end
end

