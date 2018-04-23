%% script to process reaction time data 
%% determine filename and input data
    inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Han_20180423_training\';
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
    % need to figure out exactly how I will get the true go cue times 
    % problem is that there is a delay for the stim times (for sure) and
    % possibly one for the force onset (not 100% sure)
    
    % stim can come from sync line (ainp16 typically)
    % force can come from the motor control line
    
    % the main problem is feeding them into a trial data structure
    
%% move to trial data 
    params = [];
    params.event_list = {'bumpTime'};
    td_all = parseFileByTrial(cds,params);
    
%% get movement trials
    td_move = td_all(~isnan([td_all.idx_bumpTime]));
    
%% get movement onset
    params.peak_idx_offset = 20;
    params.start_idx_offset = 15;
    params.min_ds = 1.0;
    params.which_field = 'vel';
    params.field_idx = 1;
    td_move = getMoveOnsetAndPeak(td_move,params);
    td_move = getSpeed(td_move);
%% plot a set of trials from go cue to end of trial.
    figure()
    hold on
    presample = 30;
    for t = 4:5%numel(td_move)
        yData = (td_move(t).acc(td_move(t).idx_goCueTime-presample:end,1));
        xData = ((1:numel(yData))-presample-1)*td_move(t).bin_size;
        plot(xData,yData,'k','linewidth',1.5);
        hold on
        plot((td_move(t).idx_movement_on-td_move(t).idx_goCueTime)*td_move(t).bin_size,...
            td_move(t).acc(td_move(t).idx_movement_on),...
            'r.','markersize',16)
    end
    ax = gca;
    ax.XLim = [-presample*td_move(1).bin_size,0.8];
    xlabel('Time after go cue presentation (s)')
    ylabel('Acc-x (cm/s^2)')
    ax.FontSize = 16;
    formatForLee(gcf)
%% bin the move onset 
    bE = 10:1:30;
    move_idx = [td_move.idx_movement_on] - [td_move.idx_goCueTime];
    [bC,bE] = histcounts(move_idx,bE);
    figure
    bE = bE*td_move(1).bin_size;
    bar(bE(1:end-1)+mode(diff(bE)),bC);
    
    ax = gca;
    ax.FontSize = 16;
    xlabel('RT (ms)')
    ylabel('Number of trials')
    formatForLee(gcf)


