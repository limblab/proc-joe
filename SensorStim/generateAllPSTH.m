function [] = generateAllPSTH(filepath, filename, varargin)
% to do -- deal with varargin so that the plots can be customized

%% Generates PSTH histograms for every neuron in the file one at a time
useZeroAsStart = 1;
neuronNumberStart = 1;
% load cds from file
funcFolder = pwd;

% if cds exists, load that, otherwise make new one and save cds
if(exist([filepath filename],'file') > 0)
    load([filepath filename]);
else
    labnum = 6;
    monkey = 'monkeyHan';
    ranBy = 'ranByRaeed';
    array = 'arrayLeftS1Area2';
    task = 'taskRW';
    cds = commonDataStructure();
    cds.file2cds([filepath filename], labnum, monkey, task, ranBy, array)
    save([filepath filename '_cds'],'cds');
end

cd(funcFolder);

% for each neuron, run generatePESTH

GTOstim = 0;
if(size(cds.analog,1) == 0)
    GTOstim = 1;
end
[stimState,] = determineStimTiming(cds, GTOstim, 0);
[sequenceTimes, eventTimes] = getSequenceTimes(cds, stimState,GTOstim,0);

for i = neuronNumberStart:size(cds.units,2)
    if(cds.units(i).ID ~= 0 && cds.units(i).ID ~=255) % unsorted stuff
        if(useZeroAsStart)
            generatePESTH(cds, i, 'zeroEvent', 'start','sequenceTimes',sequenceTimes,'eventTimes',eventTimes);
        else
            generatePESTH(cds, i, 'zeroEvent', 'end','sequenceTimes',sequenceTimes,'eventTimes',eventTimes);
        end
        title(['NN:' num2str(i) ' CH' num2str(cds.units(i).chan) ' ID' num2str(cds.units(i).ID)]);
        pause(3);
        close all;
    end 
end

end