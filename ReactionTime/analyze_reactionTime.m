%% script to process reaction time data 
%% determine filename and input data
    inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Han_20180426_training\';
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

    fileNumber = 1;
    cds = commonDataStructure();
    cds.file2cds([inputData.folderpath fileList(fileNumber).name],inputData.task,inputData.ranBy,...
        inputData.monkey,inputData.labnum,inputData.array1,inputData.mapFileName,'recoverPreSync');
    cd(pwd);
    % tried trial data, ran into potential issues with binning the data and
    % feeding in the true cue times. Decided to write own code
%% get cue times and information about the cue
% get movement onset, reaction time, and kinematics from just before go cue
% to end of trial for each trial
    opts.START_VAR = 'tgtOnTime';
    opts.MOVE_START_OFFSET = floor((mode([cds.trials.bumpHoldPeriod]) + 2*mode([cds.trials.bumpRisePeriod]))/mode(diff(cds.kin.t)));
    
    opts.METHOD = 'peakAcceleration';
    opts.THRESH_MULT = 0.5;
    [cueInfo] = getCueInformation(cds,opts);
    [reachData] = getReachData(cds,cueInfo,opts);

%% plot a set of reaches aligned to go cue with reaction time markers
    opts.PLOT_VAR = 'y';
    opts.NUM_PLOT = 4;
    plotReaches(reachData,cueInfo,opts);
    
%% plot a histogram showing the reaction times to each cue
    opts.BUMP_MAGS = [];
    opts.MIN_BIN = 0.1;
    opts.MAX_BIN = 0.4;
    opts.BIN_SIZE = 0.01;
    [outputData] = plotReactionTimeHistogram(reachData,cueInfo,opts);
    
%% get data related to the % success (false starts, reaches, etc.) for each cue type

    [behaviorData] = getReactionTimeSuccessRates(cds,cueInfo,opts);
