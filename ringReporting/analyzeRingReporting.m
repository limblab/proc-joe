%% set file name and load file into cds

folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\RingReporting\Chips_20180123\';
mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';

inputData.array = 'arrayLeftS1';
inputData.monkey = 'monkeyChips';
inputData.ranBy = 'ranByJoe';
inputData.lab = 6;
inputData.mapFile = strcat('mapFile',mapFileName);
inputData.task = 'taskRR';

pwd=cd;
cd(folderpath)
fileList = dir('*nev*');
cds = commonDataStructure();
cds.file2cds(strcat(folderpath,fileList(1).name),inputData.array,inputData.monkey,inputData.ranBy,...
    inputData.lab,inputData.mapFile,inputData.task,'recoverPreSync');
cd(pwd);

%%
opts.NUM_BINS_DIR = 8;
opts.MAKE_FIGURES = 1;
opts.PLOT_POLAR = 1;
opts.MAX_TRIALS_PLOT = 2000;
opts.CIRCLE_RADIUS = 7;
opts.CIRCLE_DEPTH = 2;

opts.DISTRIBUTION_BIN_SIZE = 5;
opts.BUMP_MAGS = [0.75,1.0,1.25];

behaviorData = processBehaviorRingReporting(cds,opts);

