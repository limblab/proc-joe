%% plots channels stimulated in an experiment

folderpath = 'R:\data\Han_13B1\Raw\bumpstim\20170614\proccessed\';
mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left PMd SN 6251-001460.cmp';
stimElectrodes = [35,82,57,13];
nn=1;

plotArrayMap(cds,nn,mapFileName,'numRows',10,'numCols',10,...
    'stimElectrode',stimElectrodes,'stimElectrodeColor','k','stimElectrodeLabel','string',...
    'recordingElectrode',[],'recordingElectrodeColor','k','recordingElectrodeLabel','string')
