function [times] = getTimeSinceLastSpike(spikeOccurred, pulses, currPulse)
% gets the time since the last spike has occurred for each neuron
numNeurons = size(spikeOccurred,1);

times = zeros(numNeurons,1); % initialize as zero

for idxNeuron = 1:numNeurons
    if(sum(spikeOccurred(idxNeuron,:)) > 0) % spike has occurred
        spikeIdx = 1;
        for idx = 1:1:currPulse-1
            if(spikeOccurred(idxNeuron,idx) == 1)
                spikeIdx = idx;
            end
        end
        times(idxNeuron) = pulses(currPulse,1)-pulses(spikeIdx,1);
    end
end


end