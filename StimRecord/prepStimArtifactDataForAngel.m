%% prep stim artifact data to give to Angel's group

% give them raw 30kHz data, with stim onset idx, and spike idx for easy
% data manipulation

% raw 30kHz data is in original ns5. stim onset idx is in stim info, spike
% idx is in spikes_extracted.nev

foldername = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\StimRecData\Han\Han_20190301_dblPulse\chan21stim\';
prefix = 'Han_chan21dukeProjBox_dplExp_chan21stim_A1-50_A2-50_PW1-200_PW2-200_interpulse53_pol0_002';


input_data.task='taskCObump';
input_data.ranBy='ranByJoseph'; 
input_data.array1='arrayLeftS1'; 
input_data.monkey='monkeyHan';
input_data.labnum = 6;
input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

%% load in base ns5
raw_data = openNSx([foldername,prefix,'.ns5'],'uV');

% extract raw ainp15 data
ainp15_data.data = raw_data.Data(end-1,:);
ainp15_data.t = (0:numel(ainp15_data.data)-1)/30000;
clear raw_data;

%% load in stim on times
load([foldername,prefix,'_outputData.mat']);

stimOn.idx = outputData.stimInfo.stimOn;
stimOn.ts = stimOn.idx/30000;

%% load in cds to get spike times
cds = commonDataStructure();
cds.file2cds([foldername,prefix,'_spikesExtracted'],input_data.task,input_data.ranBy,...
            input_data.monkey,input_data.labnum,input_data.array1,input_data.mapFileName); % DO NOT USE RECOVER PRE SYNC, currently this shifts the units and the analog signal differently

%% get unit ts
ts = cds.units(find([cds.units.ID] == 1)).spikes.ts;

%% package output_data
output_data.spike_times = ts;
output_data.artifact_data = ainp15_data.data;
output_data.artifact_time = ainp15_data.t;
output_data.stim_on_time = stimOn.ts;