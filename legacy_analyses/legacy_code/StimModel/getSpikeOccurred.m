function [spikeOccurred] = getSpikeOccurred(current,distances,thresholds,probabilities)
% determines if a spike occured (0 or 1) based on the current being greater
% than the threshold

% currentAdjusted = current./distances.^2;
% spikeOccurred = currentAdjusted > thresholds;

spikeOccurred = rand(size(probabilities,1),1) < probabilities(:);

end