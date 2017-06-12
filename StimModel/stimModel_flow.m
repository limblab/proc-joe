% This script runs the model

%% 1. initialize set of neurons and pulses
% a. generate list of neurons based on distances and number provided. This is
% randomely done
% b. set number of pulses (numPulses), magnitude, pulse width, and timing of
% each pulse
numNeurons = 100;
minDistance = 1; % receives I/distance^2 current
maxDistance = 3;

numPulses = 10;
magPulse = 50; % uA,10-40uA makes sense?
pulseWidth = 0.2; % ms, assume monophasic, or that only the total PW matters
timeBetween = 3; %ms

neurons = generateNeurons(numNeurons, minDistance, maxDistance);  
pulses = generatePulses(numPulses, magPulse, pulseWidth, timeBetween);

%% 2. Initialize a bunch of variables (will want to fit data later)
rheobase = 3.71; % uA
chronaxie = 0.43; % ms
tauRefractory = 112; % ms
s = 0.25; % for standard deviation in stimulation current amplitude
gamma = 2.32;
tauWindow = 40; % ms
gainThreshold = 0.4;
timeAbsolute = 1; % ms

constants = [rheobase,chronaxie,tauRefractory,s,gamma,tauWindow,gainThreshold,timeAbsolute];
%% 3. compute probability of spiking for each neuron after each pulse and 
% bump direction

[bump, spikeInfo] = computeBump(neurons,pulses,constants);


%% 4. Generate Tucker's figures by adding the stim "bump" to a physical bump at 
% different angles
