%% set initial parameters

    input_data.folderpath = 'D:\Lab\Data\BumpDirection\Han_20200819_stim_detect\';

%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\6250-002085\SN 6250-002085.cmp';
    
    input_data.date = '20190515';

    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyHan';
    input_data.ranBy = 'ranByJoe';
    input_data.lab = 6;
    input_data.mapFile = strcat('mapFile',mapFileName);
    input_data.task = 'taskBD';

    pwd=cd;
    input_data.fileNum = 1;
  
    input_data.center_y = -33;
    input_data.center_x = 3;
    input_data.num_bins = 8;
    
%% make cds
    cd(input_data.folderpath)

    file_name = dir('*nev*');
    
    params.event_list = {'goCueTime';'bumpTime';'bumpDir';'bumpMagnitude';'stimTime';'stimCode';'numTargets';'isTrainingTrial';'correctAngle'};
    params.trial_results = {'R','A','F'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [0,255];

    move_onset_params.pre_move_thresh = 1000;
    move_onset_params.min_s = 3;
    move_onset_params.start_idx_offset = -10;
    move_onset_params.max_rt_offset = 400;
    
    td_all = [];
    for i = 3%:numel(file_name)
        cds = commonDataStructure();
        cds.file2cds(strcat(input_data.folderpath,file_name(i).name),input_data.array,input_data.monkey,input_data.ranBy,...
            input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
       
        td_temp = parseFileByTrial(cds,params);
%         td_temp = stripSpikeSorting(td_temp);
        td_temp = getSpeed(td_temp);
        td_temp = getNorm(td_temp,'vel');
        td_temp = getMoveOnsetAndPeak(td_temp);
        
        td_all = [td_all,td_temp];
    end
    
    td_all = removeBadTrials(td_all);
    if(td_all(1).bin_size < 0.02)
        % set it to 20ms
        td_all = binTD(td_all,ceil(0.02/td_all(1).bin_size));
    end

%% Fix up trial data info since the kramer task doesn't output some useful things

% remove training trials and remove aborted trials where a go cue was not provided
td_use = td_all([td_all.isTrainingTrial] == 0);
td_use = td_use(~([td_use.result] == 'A' & isnan([td_use.idx_goCueTime])));
 

% get target dirs
tgt_dirs = (unique([td_use.correctAngle]));

for i_trial = 1:numel(td_use)
    % adjust bump direction to reflect true bump direction -- it comes out
    % relative to target axis
    td_use(i_trial).bumpDir = td_use(i_trial).bumpDir + td_use(i_trial).target_direction*180/pi;
    % adjust bump directions to be [0,360)
    while(td_use(i_trial).bumpDir >= 360)
        td_use(i_trial).bumpDir = td_use(i_trial).bumpDir - 360;
    end
    % set 'is_catch_trial'
    td_use(i_trial).is_catch_trial = td_use(i_trial).bumpMagnitude == 0;
    
    % get chosen target based on position at end of trial
    pos_end = td_use(i_trial).pos(td_use(i_trial).idx_endTime,:);
    pos_start = td_use(i_trial).pos(td_use(i_trial).idx_tgtOnTime,:);
    move_vec = pos_end-pos_start;
    move_angle = atan2(move_vec(2),move_vec(1));
    [~,chosen_target_idx] = min(abs(circ_dist(move_angle,tgt_dirs)));
    td_use(i_trial).chosen_target = tgt_dirs(chosen_target_idx);
end
   

%% get confusion matrix for each bump direction. Get % correct on catch trials

conf_matrix = nan(numel(tgt_dirs),numel(tgt_dirs));

% get confusion matrix
for i_pres = 1:numel(tgt_dirs)
    for i_chose = 1:numel(tgt_dirs)
         percent_presented_chosen = [td_use([td_use.correctAngle] == tgt_dirs(i_pres) & ~[td_use.is_catch_trial]).chosen_target] == tgt_dirs(i_chose);
         conf_matrix(i_pres,i_chose) = sum(percent_presented_chosen)/numel(percent_presented_chosen);
    end
end

% get percent correct on catch trials
td_catch = td_use([td_use.is_catch_trial] == 1);
percent_correct_catch = sum([td_catch.result] == 'R')/numel(td_catch);


figure();
imagesc(conf_matrix)
disp(percent_correct_catch)
    
    
%% get stim detect prob for each stim code, bumps, and catch rate
    td_stim = td_all(~isnan([td_all.stimCode]));
    td_bump = td_all([td_all.bumpMagnitude] > 0 & isnan([td_all.stimCode]));
    td_catch = td_all(isnan([td_all.stimCode]) & [td_all.bumpMagnitude] == 0 & ~isnan([td_all.idx_tgtOnTime]));

    stim_codes = unique([td_stim.stimCode]);
    percent_stim_detect = zeros(numel(stim_codes),1);
    
    for i_code = 1:numel(stim_codes)
        td_code = td_stim([td_stim.stimCode] == stim_codes(i_code));
        percent_stim_detect(i_code) = sum([td_code.result] == 'R')/numel(td_code);
    end
    
    percent_catch = sum([td_catch.result]=='R')/numel(td_catch);


