%% determine filename and input data
    input_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\CObump\Joe_20190321_moveBumps\';
%     inputData.folderpath = 'D:\Lab\Data\ReactionTime\Han_20180427_training\';
%     input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\RetiredMonkeys\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';

    input_data.task='taskCObump';
    input_data.ranBy='ranByJoseph'; 
    input_data.array1='arrayLeftS1'; 
    input_data.monkey='monkeyNone';
    input_data.labnum = 6;
    
    pwd=cd;
    cd(input_data.folderpath)
    fileList = dir('*.nev*');
    cd(pwd)
    
%% load cds, convert to td, compute PDs for all units, determine if units are well tuned
    cds = commonDataStructure();
    cds.file2cds([input_data.folderpath fileList(1).name],input_data.task,input_data.ranBy,...
        input_data.monkey,input_data.labnum,input_data.array1,input_data.mapFileName,'recoverPreSync');
    cd(pwd);
    
% convert into td
    params.event_list = {'goCueTime';'tgtDir';'bumpDir';'bumpTime'};
    params.trial_results = {'R'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [0:255];
    td_all = parseFileByTrial(cds,params);
%     td_all = removeBadTrials(td_all);
    td_all = getMoveOnsetAndPeak(td_all);
%     td_all = removeBadTrials(td_all);
%     [td_all] = removeBadNeurons(td_all,params);
    
%     % compute PD for all units
%     params = [];
%     params.array = 'LeftS1';
%     params.window = {'idx_movement_on',0;'idx_movement_on',12};
%     clear tcs
%     [tcs,confBounds,fr,covar] = getNeuronTuning(td_all,'regress',params);
%     % tcs = mean, modulation depth, pd (b0 + b1*cos(x-b2))
%     tcs = tcs(:,1:3);
%     
%     % determine if units are well tuned (conf bounds less than +- 45 deg)
%     pd_idx = 3;
%     conf_diff = angleDiff(confBounds{pd_idx}(:,1),tcs(:,pd_idx),1,0);
%     well_tuned = abs(conf_diff) < pi/2;

%% seperate trials
td_perp_bump = td_all(~isnan([td_all.idx_bumpTime]) & [td_all.bumpDir] == 90);
td_parallel_bump = td_all(~isnan([td_all.idx_bumpTime]) & ([td_all.bumpDir] == 0 | [td_all.bumpDir] == 180));

td_noBump = td_all(isnan([td_all.idx_bumpTime]) & [td_all.tgtDir] == 0);

%% plot reaches to targets
f=figure();
f.Name = 'No bump';
color_to_use = [linspace(0,0.75,10)',linspace(0,0.75,10)',linspace(0,0.75,10)'];
counter = 1;
for tr = 1:15%numel(td_noBump)
    window = [td_noBump(tr).idx_goCueTime,td_noBump(tr).idx_endTime+20];
    reach_data = td_noBump(tr).pos(window(1):window(2),:);
    if(any(reach_data(:,2) > -31) || any(reach_data(:,2) < -32.5))
    else
        plot(reach_data(:,1),reach_data(:,2),'color',color_to_use(counter,:),'linewidth',2);
        hold on
        counter = counter+1;
    end
end

% f=figure();
% f.Name = 'Perpendicular_bumps';
num_trials = 40;
color_to_use = [1+zeros(num_trials,1),linspace(0,0.5,num_trials)',linspace(0,0.5,num_trials)'];
    
for tr = 1:numel(td_perp_bump)
    figure();
    window = [td_perp_bump(tr).idx_goCueTime,td_perp_bump(tr).idx_endTime+5];
    reach_data = td_perp_bump(tr).pos(window(1):window(2),:);
    
    plot(reach_data(:,1),reach_data(:,2),'color',color_to_use(tr,:),'linewidth',2);
    hold on
    plot(td_perp_bump(tr).pos(td_perp_bump(tr).idx_bumpTime,1),td_perp_bump(tr).pos(td_perp_bump(tr).idx_bumpTime,2),...
        '.','color',getColorFromList(1,2),'markersize',24)
    plot(td_perp_bump(tr).pos(td_perp_bump(tr).idx_bumpTime+20,1),td_perp_bump(tr).pos(td_perp_bump(tr).idx_bumpTime+20,2),...
        '.','color',getColorFromList(1,2),'markersize',24)
    plot(td_perp_bump(tr).pos(td_perp_bump(tr).idx_bumpTime+40,1),td_perp_bump(tr).pos(td_perp_bump(tr).idx_bumpTime+40,2),...
        '.','color',getColorFromList(1,3),'markersize',24)
    plot(td_perp_bump(tr).pos(td_perp_bump(tr).idx_bumpTime+55,1),td_perp_bump(tr).pos(td_perp_bump(tr).idx_bumpTime+55,2),...
        '.','color',getColorFromList(1,3),'markersize',24)
end

% f=figure();
% f.Name = 'Parallel_bumps';
% for tr = 1:numel(td_parallel_bump)
%     window = [td_parallel_bump(tr).idx_goCueTime,td_parallel_bump(tr).idx_endTime];
%     reach_data = td_parallel_bump(tr).pos(window(1):window(2),:);
%     plot(reach_data(:,1),reach_data(:,2),'r');
%     hold on
% end

%% distribution of theta
no_bump_pos = [];
bump_pos = [];
window = [40,55];
theta_no_bump = [];
for tr = 1:numel(td_noBump)
    no_bump_pos(tr,:,:) = td_noBump(tr).pos(td_noBump(tr).idx_goCueTime:td_noBump(tr).idx_goCueTime+50,:);
end

for tr = 1:numel(td_perp_bump)
    bump_pos(tr,:,:) = td_perp_bump(tr).pos(td_perp_bump(tr).idx_bumpTime+window(1):td_perp_bump(tr).idx_bumpTime+window(2),:);
end

mean_no_bump = squeeze(mean(no_bump_pos,1));
mean_no_bump_vec = mean_no_bump(end,:) - mean_no_bump(1,:);
mean_no_bump_theta = atan2(mean_no_bump_vec(2),mean_no_bump_vec(1));

no_bump_thetas = 180/pi*(atan2(no_bump_pos(:,end,2)-no_bump_pos(:,1,2),no_bump_pos(:,end,1)-no_bump_pos(:,1,1)) - mean_no_bump_theta);
bump_thetas = 180/pi*(atan2(bump_pos(:,end,2)-bump_pos(:,1,2),bump_pos(:,end,1)-bump_pos(:,1,1))-mean_no_bump_theta);

figure
bin_edges = -180:10:180;
no_bump_count = histcounts(no_bump_thetas,bin_edges);
bump_count = histcounts(bump_thetas,bin_edges);

histogram(no_bump_thetas,bin_edges,'FaceColor','k')
hold on
histogram(bump_thetas,bin_edges,'FaceColor','r')
%% velocity profile

f=figure();
% f.Name = 'No bump';
for tr = 1:5%:numel(td_noBump)
    window = [td_noBump(tr).idx_goCueTime,td_noBump(tr).idx_endTime+8];
    reach_data = td_noBump(tr).acc(window(1):window(2),:);
    plot(reach_data(:,2),'k');
    hold on
end

% f=figure();
% f.Name = 'Perpendicular_bumps';
for tr = 1:5%numel(td_perp_bump)
    window = [td_perp_bump(tr).idx_goCueTime,td_perp_bump(tr).idx_endTime+8];
    reach_data = td_perp_bump(tr).acc(window(1):window(2),:);
    plot(reach_data(:,2),'b');
    hold on
    plot(td_perp_bump(tr).pos(td_perp_bump(tr).idx_bumpTime,1),td_perp_bump(tr).pos(td_perp_bump(tr).idx_bumpTime,2),...
        'r.','markersize',20)
    plot(td_perp_bump(tr).pos(td_perp_bump(tr).idx_bumpTime+20,1),td_perp_bump(tr).pos(td_perp_bump(tr).idx_bumpTime+20,2),...
        'r.','markersize',20)
end

% f=figure();
% % f.Name = 'Parallel_bumps';
% for tr = 1:10%numel(td_parallel_bump)
%     if(td_parallel_bump(tr).bumpDir ~= td_parallel_bump(tr).tgtDir)
%         window = [td_parallel_bump(tr).idx_goCueTime,td_parallel_bump(tr).idx_endTime];
%         reach_data = td_parallel_bump(tr).vel(window(1):window(2),:);
%         plot(sum(reach_data.^2,2),'b');
%         hold on
%     end
% end


%% change point analysis...
td_noBump = td_all(isnan([td_all.idx_bumpTime]));
trial = 2;
% plot(td_noBump(trial).pos(td_noBump(trial).idx_goCueTime:td_noBump(trial).idx_endTime,1),...
%     td_noBump(trial).pos(td_noBump(trial).idx_goCueTime:td_noBump(trial).idx_endTime,2))

reach_data = td_noBump(trial).pos(td_noBump(trial).idx_goCueTime:td_noBump(trial).idx_endTime,2);
plot(reach_data)
ipt = findchangepts(reach_data)