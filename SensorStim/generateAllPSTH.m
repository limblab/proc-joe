function [neuronRet] = generateAllPSTH(filepath, filename, varargin)
% to do -- deal with varargin so that the plots can be customized

%% Generates PSTH histograms for every neuron in the file one at a time
useZeroAsStart = 1;
neuronNumberStart = 1;
% load cds from file
funcFolder = pwd;

noPlots = 0;
useRate = 1;
ret = 0;
binSize = 0.05;
spindleStim = 0; % GTO stim is 0, spindle stim is 1
neuronRet = [];
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'noPlots' % so nothing is plotted
            noPlots = varargin{i+1};
        case 'useRate' % returns things in terms of rates
            useRate = varargin{i+1};
        case 'return' % returns a set of neurons that matter
            ret = varargin{i+1};
        case 'binSize'
            binSize = varargin{i+1};
        case 'spindleStim'
            spindleStim = varargin{i+1};
    end
    
end

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
if((~exist('sequenceTimes') || ~exist('eventTimes')) || ~exist('stimState'))
    [stimState,] = determineStimTiming(cds, GTOstim, 0);
    [sequenceTimes, eventTimes] = getSequenceTimes(cds, stimState,GTOstim,0);
    save([filepath filename],'cds','sequenceTimes','eventTimes','stimState');
end

for i = neuronNumberStart:size(cds.units,2)
    if(cds.units(i).ID ~= 0 && cds.units(i).ID ~=255) % unsorted stuff
        if(useZeroAsStart)
            [binEdges,binCounts]=generatePESTH(cds, i, 'zeroEvent', 'start','sequenceTimes',sequenceTimes,'eventTimes',eventTimes,...
                'stimState',stimState,'noPlots',noPlots,'useRate',useRate, 'return', ret,'binSize',binSize);
        else
            [binEdges,binCounts]=generatePESTH(cds, i, 'zeroEvent', 'end','sequenceTimes',sequenceTimes,'eventTimes',eventTimes,...
                'stimState',stimState,'noPlots',noPlots,'useRate',useRate, 'return', ret,'binSize',binSize);
        end
        if(~noPlots)
            title(['NN:' num2str(i) ' CH' num2str(cds.units(i).chan) ' ID' num2str(cds.units(i).ID)]);
            pause(3);
            close all;
        end
        if(~spindleStim && ret && keepNeuronGTOstim(cds, i, binEdges, binCounts, eventTimes(1),useRate)==1) % determine if keeping
            neuronRet = [neuronRet; i];
        end
        if(spindleStim && ret && keepNeuronSpindleStim(cds, i, binEdges, binCounts, eventTimes,useRate,stimState)==1) % determine if keeping
            neuronRet = [neuronRet; i];
        end
    end 
end

end