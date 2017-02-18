%% this script is similar to generateAllPSTH -- except it rejects neurons that
% do not meet some threshold. The goal is to speed up the process of
% looking for neurons by autorejecting some :D

% load cds from file
funcFolder = pwd;

filepath = 'D:\Lab\Data\SensorStim\Han_20170209\';
filename = 'Han_20170209_GTOstim_FCU_1mA_area2_001';

% if cds exists, load that, otherwise make new one and save cds
if(exist([filepath filename '_cds.mat'],'file') > 0)
    load([filepath filename '_cds.mat']);
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

% for each neuron, run generatePESTH to get bincounts and binedges,
% suppress output

GTOstim = 0;
if(size(cds.analog,1) == 0)
    GTOstim = 1;
end
[stimState,] = determineStimTiming(cds, GTOstim, 0);
[sequenceTimes, eventTimes] = getSequenceTimes(cds, stimState,GTOstim,0);
out_mat = [];
for i = 1:size(cds.units,2)
    if(cds.units(i).ID ~= 0 && cds.units(i).ID ~=255) % unsorted stuff
        if(useZeroAsStart)
            [bE,bC]=generatePESTH(cds, i, 'zeroEvent', 'start','sequenceTimes',sequenceTimes,...
                'noPlots',1,'eventTimes',eventTimes);
        else
            [bE,bC]=generatePESTH(cds, i, 'zeroEvent', 'end','sequenceTimes',sequenceTimes,... 
                'noPlot',1,'eventTimes',eventTimes);
        end
        
        if(keepNeuronGTOstim(cds,i,bE,bC)==1)
            % store the neuron in an output matrix
            out_mat = [out_mat i];
        end
        
    end 
end