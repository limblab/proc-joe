%% script to process reaction time data 
%% determine filename and input data
    inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Han_20180418_training\';
    inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

    inputData.task='taskRT';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayLeftS1'; 
    inputData.monkey='monkeyHan';
    inputData.labnum = 6;
    
    pwd=cd;
    cd(inputData.folderpath)
    fileList = dir('*.nev*');

%% load in cds

    cds = commonDataStructure();
    cds.file2cds([inputData.folderpath fileList(1).name],inputData.task,inputData.ranBy,...
        inputData.monkey,inputData.labnum,inputData.array1,inputData.mapFileName,'recoverPreSync');

    % need to figure out exactly how I will get the true go cue times 
    % problem is that there is a delay for the stim times (for sure) and
    % possibly one for the force onset (not 100% sure)
    
    % stim can come from sync line (ainp16 typically)
    % force can come from the motor control line
    
%% move to trial data 
    params = [];
%     params.event_names = {''};
    td_all = parseFileByTrial(cds,params);