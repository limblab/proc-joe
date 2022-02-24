function [ bump, spikeInfo ] = computeBump(neurons, pulses, constants )

% extract constants
rheobase = constants(1); % uA
chronaxie = constants(2); % ms
tauRefractory = constants(3); % ms
s = constants(4); % for standard deviation in stimulation current amplitude
gamma = constants(5);
tauWindow = constants(6); % ms
gainThreshold = constants(7);
timeAbsolute = constants(8); % ms

numNeurons = size(neurons,1);
%% Compute probability of spiking for each neuron after each pulse
% look at getSpikeOccurred -- not sure if doing correctly

spikeProbability = [];
spikeOccurred = [];
thresholds = [];
for pulse = 1:size(pulses,1) % for each pulse
    if(pulse == 1)
        time = zeros(numNeurons,1); % first pulse at time = 0, no spike occurs before this
    else
        % if a spike has occurred, time is the time since that spike.
        % Otherwise, time is equal to 0 (no spike has occurred)
        time = getTimeSinceLastSpike(spikeOccurred, pulses, pulse);
    end
    
    % get threshold based on time of last spike
    for idxNeuron = 1:numNeurons
        thresholds(idxNeuron,1) = getThreshold(time(idxNeuron),timeAbsolute,tauRefractory, gamma, pulses(pulse,3), chronaxie, rheobase);
    end
    
    spikeProbability(:,pulse) = probSpike(pulses(pulse,2),neurons(:,1),thresholds, gainThreshold, s);
    spikeOccurred(:,pulse) = getSpikeOccurred(pulses(pulse,2),neurons(:,1),thresholds, spikeProbability(:,pulse));
end

%% Based on spike data, time window thing, and neuron directions, compute a vector "bump"
timePoints = pulses(:,1);
timeWindow = exp(-timePoints/tauWindow);
spikePower = spikeProbability*timeWindow;
spikeDirections = [spikePower.*neurons(:,2), spikePower.*neurons(:,3)];
bump = sum(spikeDirections);

spikeInfo.spikeProbability = spikeProbability;
spikeInfo.spikeOccurred = spikeOccurred;
spikeInfo.spikePower = spikePower;
spikeInfo.spikeDirections = spikeDirections;
end

