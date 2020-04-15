%% Reza's emg result
%% determine filename and input data
    inputData.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\CObump\';

    inputData.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    
    
    inputData.task='taskCObump';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayLeftS1'; 
    inputData.monkey='monkeyHan';
    inputData.labnum = 6;
    
    inputData.center = [-2,2;-36,-34]; % x bounds, y bounds
    
    pwd=cd;
    cd(inputData.folderpath)
    fileList = dir('*.nev*');
    cd(pwd)

%% make cds
    cds = commonDataStructure();
    cds.file2cds([inputData.folderpath fileList(1).name],inputData.task,inputData.ranBy,...
        inputData.monkey,inputData.labnum,inputData.array1,inputData.mapFileName);
    cd(pwd);

    
%% convert to td

    params.event_list = {'tgtOnTime';'bumpTime';'bumpMagnitude';'bumpDir';'tgtDir'};
    params.trial_results = {'R','F'};
    params.all_points = 1;
    
    params.include_ts = 0;
    params.exclude_units = [0,255];
    td_all = parseFileByTrial(cds,params);
    
%     td_no_bump = td_all(isnan([td_all.bumpDir]));
%     td_no_bump = td_no_bump([td_no_bump.result] == 'R');

%% concatenate TD
    td_cat = [];
    td_cat.pos = []; td_cat.vel = [];
    td_cat.acc = []; td_cat.force = []; td_cat.emg = [];
    for tr = 1:numel(td_all)
        td_cat.pos = [td_cat.pos;td_all(tr).pos];
        td_cat.vel = [td_cat.vel;td_all(tr).vel];
        td_cat.acc = [td_cat.acc;td_all(tr).acc];
        td_cat.force = [td_cat.force;td_all(tr).force];
        td_cat.emg = [td_cat.emg;td_all(tr).emg];
    end
    
    td_cat.target_direction = [td_all.target_direction];
    td_cat.trial_id = [td_all.trial_id];
    td_cat.bin_size = [td_all.bin_size];
    td_cat.bumpDir = [td_all.bumpDir];
    td_cat.bumpMagnitude = [td_all.bumpMagnitude];
    td_cat.tgtDir = [td_all.tgtDir];
    td_cat.idx_startTime = [td_all.idx_startTime];
    td_cat.idx_bumpTime = [td_all.idx_bumpTime];
    td_cat.idx_tgtOnTime = [td_all.idx_tgtOnTime];
    td_cat.idx_goCueTime = [td_all.idx_goCueTime];
    td_cat.idx_endTime = [td_all.idx_endTime];
    td_cat.emg_names = td_all.emg_names;
    
%% determine when animal is still
% call him still if speed is smaller than some threshold
    is_still_speed_thresh = 1;

    td_cat.speed = sqrt(td_cat.vel(:,1).^2 + td_cat.vel(:,2).^2);
    td_cat.is_still = td_cat.speed < is_still_speed_thresh;
    td_cat.is_center = td_cat.pos(:,1) > inputData.center(1,1) & td_cat.pos(:,1) < inputData.center(1,2) & ...
        td_cat.pos(:,2) > inputData.center(2,1) & td_cat.pos(:,2) < inputData.center(2,2);
    % plot cursor position with when still as a different color
    
    idx_start = td_cat.idx_goCueTime(1)-200;
    idx_end = td_cat.idx_endTime(1)+2000;
    still_mask = td_cat.is_still(idx_start:idx_end) & td_cat.is_center(idx_start:idx_end);
    
    
    pos_plot = td_cat.pos(idx_start:idx_end,:);
    figure();
    plot(pos_plot(still_mask == 0,1),pos_plot(still_mask == 0,2),'.')
    hold on
    plot(pos_plot(still_mask == 1,1),pos_plot(still_mask == 1,2),'r.','markersize',20)
    

%% use is_still and is_center to get reaches from center to outer to center
    min_reach_length = 6; % all reaches must have length > 6
    
    
    


    

