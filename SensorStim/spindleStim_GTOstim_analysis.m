%%%% Raster plot of all spikes for the trial

sortThreshold = 0;

%% get data into cds format
funcFolder = pwd;

filepath = 'D:\Lab\Data\SensorStim\Chips_20151123\';
% filename = 'Chips_20151123_TricepsSweep_009'; % NN: 36
% filename = 'Chips_20151123_WristExtVibe_007'; % NN: 3, 93, 4
% filename = 'Chips_20151123_WristFlexVibe_006'; % NN: 21

filename = 'Chips_20151123_GTOStim_03mA_noartefactrejection_003';

% filepath = 'D:\Lab\Data\SensorStim\Han_20170106\';
% filename = 'Han_20170106_SpindleStim_FCR_area2EMG_003'; % 156
% filename = 'Han_20170106_SpindleStim_FCRdistal_area2EMG_004'; % 59, 109

GTOstim = ~isempty(strfind(lower(filename),'gtostim'));

timeAfterGTOStim = 0.5; % in seconds

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

%% raster: all units against time sorted by firing frequency ratio if that ratio is greater than sortThreshold
[f,sortFFR] = plotRaster(cds, sortThreshold, GTOstim, timeAfterGTOStim);
