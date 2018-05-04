%% script to process reaction time data 
%% determine filename and input data
    inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Han_20180501\';
    inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

    inputData.task='taskRT';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayLeftS1'; 
    inputData.monkey='monkeyHan';
    inputData.labnum = 6;
    
    pwd=cd;
    cd(inputData.folderpath)
    fileList = dir('*.nev*');

%% load in cds and extract data

    cds = commonDataStructure();
    cds.file2cds([inputData.folderpath fileList(fileNumber).name],inputData.task,inputData.ranBy,...
        inputData.monkey,inputData.labnum,inputData.array1,inputData.mapFileName,'recoverPreSync');
    cd(pwd);
    
%% plot a set of reaches aligned to go cue with reaction time markers

    
%% plot a histogram showing the reaction times to each cue

    
%% plot the mean reaction time to each cue based on magnitude

    
%% plot a psychometic curve fitting detection probability for the provided data
    
    
    
    
    