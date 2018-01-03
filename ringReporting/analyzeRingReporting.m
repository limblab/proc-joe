%% set file name and load file into cds

folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\RingReporting\Chips_20171228\';
mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';

inputData.array = 'arrayLeftS1';
inputData.monkey = 'monkeyChips';
inputData.ranBy = 'ranByJoe';
inputData.lab = 6;
inputData.mapFile = strcat('mapFile',mapFileName);
inputData.task = 'taskUNT2D';

pwd=cd;
cd(folderpath)
fileList = dir('*-s*');
cds = commonDataStructure();
cds.file2cds(strcat(folderpath,fileList(1).name),inputData.array,inputData.monkey,inputData.ranBy,...
    inputData.lab,inputData.mapFile,inputData.task);
cd(pwd);

%%
opts.NUM_BINS_DIR = 8;
opts.MAKE_FIGURES = 0;
opts.PLOT_POLAR = 0;

behaviorData = processBehaviorRingReporting(cds,opts);
