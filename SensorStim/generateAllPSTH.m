%% Generates PSTH histograms for every neuron in the file one at a time
useZeroAsStart = 1;
neuronNumberStart = 1;
% load cds from file
funcFolder = pwd;

filepath = 'D:\Lab\Data\SensorStim\Chips_20151123\';
% filename = 'Chips_20151123_GTOStim_03mA_artefactrejection_002'; % 17, 21, 36, 37, 84, 118, 136, 183, 189
% filename = 'Chips_20151123_GTOStim_03mA_noartefactrejection_003'; % 7,13,36,40,60,65,87,96,119
% filename = 'Chips_20151123_GTOStim_025mA_noartefactrejection_001'; % 67,101,193

% filepath = 'D:\Lab\Data\SensorStim\Han_20170106\';
% filename = 'Han_20170106_SpindleStim_FCR_area2EMG_003'; 

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
        pause(5);
        close all;
    end 
end