%% use this script to design biomimetic stimulation patterns


%% determine filename and input data
    input_data.folderpath = 'C:\Users\jts3256\Desktop\Han_CObump\';
%     inputData.folderpath = 'D:\Lab\Data\ReactionTime\Han_20180427_training\';
%     input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\RetiredMonkeys\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
    
    input_data.task='taskCObump';
    input_data.ranBy='ranByJoseph'; 
    input_data.array1='arrayLeftS1'; 
    input_data.monkey='monkeyHan';
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
    
%% convert into td
    params.event_list = {'goCueTime';'tgtDir';'bumpDir';'bumpTime'};
    params.trial_results = {'R'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [1,255];
    td_all = parseFileByTrial(cds,params);
    td_all = removeBadTrials(td_all);
    td_all = getMoveOnsetAndPeak(td_all);
    td_all = removeBadTrials(td_all);
    [td_all] = removeBadNeurons(td_all,params);
    
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
    
%% for each neuron, plot the mean firing rate for each unique BD (like a PD plot, but no cosine tuning)
    bump_dir = unique([td_all.bumpDir]);
    mean_fr_bump = [];
    mean_fr_baseline = [];
    chan_num = [];
    for unit = 1:10%size(td_all(1).LeftS1_unit_guide,1)
        for bd = 1:numel(bump_dir)
            trial_mask = [td_all.bumpDir] == bump_dir(bd);
            fr_data_bump = [];
            fr_data_baseline = [];
            td_temp = td_all(trial_mask);
            for t = 1:numel(td_temp)
                window_bump = td_temp(t).idx_bumpTime + [0,12];
                window_baseline = td_temp(t).idx_bumpTime - [14,2];
%                 window_baseline = [td_all(t).idx_startTime, td_all(t).idx_bumpTime-2];
                fr_data_bump(t) = sum(td_temp(t).LeftS1_spikes(window_bump(1):window_bump(2),unit))/(diff(window_bump)*td_temp(t).bin_size);
                fr_data_baseline(t) = sum(td_temp(t).LeftS1_spikes(window_baseline(1):window_baseline(2),unit))/(diff(window_baseline)*td_temp(t).bin_size);
            end
            mean_fr_bump(unit,bd) = mean(fr_data_bump);
            mean_fr_baseline(unit,bd) = mean(fr_data_baseline);
        end
        chan_num(unit,1) = td_all(1).LeftS1_unit_guide(unit,1);

        figure()
        plot(bump_dir,mean_fr_bump(unit,:) - mean(mean_fr_baseline(unit,:),2),'.','markersize',14)
    end

mean_fr = mean_fr_bump - mean(mean_fr_bump,2);

%% design stim pattern
mean_fr_norm = mean_fr;
mean_fr_norm(mean_fr < 0) = 0;
mean_fr_norm = mean_fr_norm./max(mean_fr_norm,[],2);

%% plot heatmap of stim pattern
map_filename = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
map_data = loadMapFile(map_filename);
colors = [(0:1:63)'/63 zeros(64,1) zeros(64,1)];

for bd = 1:size(mean_fr_norm,2)
    figure();
    plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5)
    hold on
    for elec = 1:size(chan_num,1)
        map_idx = find(chan_num(elec) == map_data.chan);
        freq_norm = mean_fr_norm(elec,bd);
        if(freq_norm > 0)
            color_idx = max(1,min(size(colors,1),ceil(freq_norm*size(colors,1))));

            rectangle('Position',[map_data.row(map_idx),map_data.col(map_idx),1,1],'FaceColor',colors(color_idx,:), 'EdgeColor',colors(color_idx,:));
            hold on
        end
    end
    
    f=gcf;
    set(gca,'visible','off')
%     f.Name = ['Han_20180802_BD_heatmap_',num2str(bump_dir(bd))];
%     saveFiguresLIB(f,input_data.folderpath,f.Name);
end

%% put desired axis into wave_mappings
    axis_idx = 3;
    wave_mappings = {};
    wave_mappings{1}(:,:) = [chan_num,mean_fr_norm(:,axis_idx)];
    wave_mappings{2}(:,:) = [chan_num,mean_fr_norm(:,axis_idx+4)];

    wave_mappings{1}(wave_mappings{1}(:,2)<0.15,:) = [];
    wave_mappings{2}(wave_mappings{2}(:,2)<0.15,:) = [];

%% dissimilarity metric
    % euclidean distance between electrode and nearest electrode multiplied
    % by normalized freq. 


