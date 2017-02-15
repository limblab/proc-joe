%% Generates PSTH histograms for every neuron in the file one at a time
useZeroAsStart = 1;
neuronNumberStart = 40;
% load cds from file
funcFolder = pwd;

% filepath = 'D:\Lab\Data\SensorStim\Chips_20151123\';
% filename = 'Chips_20151123_GTOStim_03mA_artefactrejection_002'; % 21, 36, 136
% filename = 'Chips_20151123_GTOStim_03mA_noartefactrejection_003'; % 7,36,119
% filename = 'Chips_20151123_GTOStim_025mA_noartefactrejection_001'; % 67,101,193

filepath = 'D:\Lab\Data\SensorStim\Han_20170209\';
filename = 'Han_20170209_GTOstim_ECR_2mA_area2_001';

GTOstim = ~isempty(strfind(lower(filename),'gtostim'));

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

% for each neuron, run generatePESTH
for i = neuronNumberStart:size(cds.units,2)
    if(cds.units(i).ID ~= 0 && cds.units(i).ID ~=255) % unsorted stuff
        if(useZeroAsStart)
            generatePESTH(cds, i, 'zeroEvent', 'start');
        else
            generatePESTH(cds, i, 'zeroEvent', 'end');
        end
        title(['NN:' num2str(i) ' CH' num2str(cds.units(i).chan) ' ID' num2str(cds.units(i).ID)]);
        pause(3);
        close all;
    end 
end