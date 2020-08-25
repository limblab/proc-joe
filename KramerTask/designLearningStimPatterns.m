%% use this script to design biomimetic stimulation patterns


%% determine filename and input data
    input_data.folderpath = 'D:\Lab\Data\BumpDirection\Han_20200820_cobump\';
%     input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    input_data.mapFileName = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

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
    
  % covert to td
    params.event_list = {'goCueTime';'bumpTime';'bumpDir';'bumpMagnitude';'stimTime';'stimCode';'numTargets';'isTrainingTrial';'correctAngle'};
    params.trial_results = {'R'};
    params.extra_time = [1,2];
    params.include_ts = 1;
    params.exclude_units = [255];
    
    params.remove_nan_idx = true;
    params.nan_idx_names = {'idx_goCueTime','idx_startTime','idx_endTime','idx_movement_on'};
    td = parseFileByTrial(cds,params);
    td = getNorm(td,{'vel'});
    td = stripSpikeSorting(td);

    params.start_idx = 'idx_goCueTime';
    params.end_idx = 'idx_endTime';
    td = getMoveOnsetAndPeak(td,params);
    td = removeBadTrials(td,params);

    % remove training trials and remove aborted trials where a go cue was not provided
%     td = td([td.isTrainingTrial] == 0);
 
%     for i_trial = 1:numel(td)
%         % adjust bump direction to reflect true bump direction -- it comes out
%         % relative to target axis
%         td(i_trial).bumpDir = td(i_trial).bumpDir + td(i_trial).target_direction*180/pi;
%         % adjust bump directions to be [0,360)
%         while(td(i_trial).bumpDir >= 360)
%             td(i_trial).bumpDir = td(i_trial).bumpDir - 360;
%         end
%     end
    
    % get trials with center hold bumps only
    td_bump_all = td([td.idx_goCueTime] > [td.idx_bumpTime] & [td.bumpMagnitude] > 0);
    
    for tr = 1:numel(td_bump_all)
        td_bump_all(tr).bumpDir = mod(td_bump_all(tr).bumpDir,360);
    end
    
    td_move_all = td(isnan([td.idx_bumpTime]) | [td.idx_goCueTime] > [td.idx_bumpTime]);
    
    chan_list_all = td(1).LeftS1_unit_guide(:,1);
%% remove channels not on chan_list from trial data

    keep_mask = ones(size(td_bump_all(1).LeftS1_unit_guide,1),1);
    for i_unit = 1:numel(keep_mask)
        keep_mask(i_unit) = sum(chan_list == td_bump_all(1).LeftS1_unit_guide(i_unit,1)) > 0;
    end

    for i_trial = 1:numel(td_bump_all)
        td_bump_all(i_trial).LeftS1_spikes = td_bump_all(i_trial).LeftS1_spikes(:,keep_mask==1);
        td_bump_all(i_trial).LeftS1_unit_guide = td_bump_all(i_trial).LeftS1_unit_guide(keep_mask == 1,:);
    end
    
    keep_mask = ones(size(td_move_all(1).LeftS1_unit_guide,1),1);
    for i_unit = 1:numel(keep_mask)
        keep_mask(i_unit) = sum(chan_list == td_move_all(1).LeftS1_unit_guide(i_unit,1)) > 0;
    end
    for i_trial = 1:numel(td_move_all)
        td_move_all(i_trial).LeftS1_spikes = td_move_all(i_trial).LeftS1_spikes(:,keep_mask==1);
        td_move_all(i_trial).LeftS1_unit_guide = td_move_all(i_trial).LeftS1_unit_guide(keep_mask == 1,:);
    end

    chan_list_all = td_bump_all(1).LeftS1_unit_guide(:,1);
 %% plot raster (or PSTH) during bumps in each direction
 
    desired_chan = 14;
    i_unit = find(td_bump_all(1).LeftS1_unit_guide(:,1) == desired_chan);
%     i_unit = 49;
    window_ms = [-250,500];
    window_idx = (window_ms/1000)/td_bump_all(1).bin_size;
        
    bump_dirs = unique([td_bump_all.bumpDir]);
%     for i_unit = 1:numel(td_bump_all(1).LeftS1_ts)
        figure('Position',[680 120 560 858]); hold on;
        for i_dir = 1:numel(bump_dirs)
            x_data = []; y_data = [];
            trial_counter = 1;
            td_bump_dir = td_bump_all([td_bump_all.bumpDir] == bump_dirs(i_dir));
            
            for i_trial = 1:numel(td_bump_dir)
                spike_data = td_bump_dir(i_trial).LeftS1_ts{i_unit} - (td_bump_dir(i_trial).idx_bumpTime - td_bump_dir(i_trial).idx_startTime)*td_bump_dir(i_trial).bin_size;
                spike_mask = spike_data > window_ms(1)/1000 & spike_data < window_ms(2)/1000;
                y_data(end+1:end+sum(spike_mask)) = i_trial;
                x_data(end+1:end+sum(spike_mask)) = spike_data(spike_mask);
            end
            subplot(8,1,i_dir)
            plot(x_data,y_data,'k.','markersize',2)
            ylim([0,i_trial])
            xlim(window_ms/1000)
        end
%     end
    
    
%% plot psth for movement around movement onset
    move_window_ms = [-500,1000];
    move_window_idx = (move_window_ms/1000)/td_bump_all(1).bin_size;
    td_move = trimTD(td_move_all,{'idx_movement_on',move_window_idx(1)},{'idx_movement_on',move_window_idx(2)});

    td_move = binTD(td_move,20);
    
    move_dirs = unique([td_move.target_direction]);
    for i_unit = 1:size(td_move(1).LeftS1_spikes,1)
        figure('Position',[680 120 560 858]);
        ax_list = [];
        for i_dir = 1:numel(move_dirs)
            trial_mask = isEqual([td_move.target_direction], move_dirs(i_dir));
            td_move_dir = td_move(trial_mask==1);
            td_move_avg = trialAverage(td_move_dir);
            idx_move = td_move_avg(1).idx_movement_on;
            
            ax_list(end+1) = subplot(8,1,i_dir);
            plot((0:size(td_move_avg(1).LeftS1_spikes,1)-1) - idx_move, td_move_avg(1).LeftS1_spikes(:,i_unit))
        end
        linkaxes(ax_list,'xy');
        xlim([-idx_move,size(td_move_avg(1).LeftS1_spikes,1)-1-idx_move]);
    end
    
%% get PD for each neuron and mean FR during bumps in all tested directions
    bump_window_ms = [-500,200];
    bump_delay_ms = 80; % delay before counting spikes post bump
    
    delay_offset = floor(bump_delay_ms/td_bump_all(1).bin_size/1000);
    bump_window_idx = (bump_window_ms/1000)/td_bump_all(1).bin_size;
    
    td_bump = trimTD(td_bump_all,{'idx_bumpTime',bump_window_idx(1)},{'idx_bumpTime',bump_window_idx(2)});
    
    move_window_ms = [-100,250];
    move_window_idx = (move_window_ms/1000)/td_bump_all(1).bin_size;
    td_move = trimTD(td_move_all,{'idx_movement_on',move_window_idx(1)},{'idx_movement_on',move_window_idx(2)});
    
    % get PDs
    pd_params = [];
    pd_params.out_signals = 'LeftS1_spikes';
    pd_params.bootForTuning = 0;
    pd_params.num_boots = 1;
    pd_params.move_corr = 'vel';
    
    pd_bump = getTDPDs(td_bump,pd_params);
    pd_move = getTDPDs(td_move,pd_params);
 %   
    % get mean FR during bump for each dir
    bump_dirs = unique([td_bump.bumpDir]);
    num_units = size(chan_list_all,1);
    mean_fr_bump = zeros(numel(bump_dirs),num_units);
    mean_fr_move = zeros(size(mean_fr_bump));
    mean_fr_baseline = zeros(size(mean_fr_bump));
    std_fr_baseline = zeros(size(mean_fr_bump));
    z_score_all = zeros(size(mean_fr_bump));
    
    for i_bump = 1:numel(bump_dirs)
        trial_mask = isEqual(round([td_bump.bumpDir]), bump_dirs(i_bump));
        td_bump_dir = td_bump(trial_mask==1);
        trial_mask = isEqual(round([td_move.target_direction]*180/pi), bump_dirs(i_bump));
        td_move_dir = td_move(trial_mask==1);
        td_bump_avg = trialAverage(td_bump_dir);
        td_move_avg = trialAverage(td_move_dir);
        idx_bump = td_bump_avg(1).idx_bumpTime;
        idx_move = td_move_avg(1).idx_movement_on;
        
        mean_fr_move(i_bump,:) = mean(td_move_avg(1).LeftS1_spikes(idx_move:end))/td_move_avg(1).bin_size;
        mean_fr_bump(i_bump,:) = mean(td_bump_avg(1).LeftS1_spikes(idx_bump+delay_offset:end,:))/td_bump_avg(1).bin_size;
        mean_fr_baseline(i_bump,:) = mean(td_bump_avg(1).LeftS1_spikes(1:idx_bump-10,:))/td_bump_avg(1).bin_size;
        fr_temp = zeros(numel(td_bump_dir),num_units);
        for i_trial = 1:numel(td_bump_dir)
            fr_temp(i_trial,:) = mean(td_bump_dir(i_trial).LeftS1_spikes(1:idx_bump-10,:))/td_bump_avg(1).bin_size;
        end
        std_fr_baseline(i_bump,:) = std(fr_temp);
    end
    std_fr_baseline(std_fr_baseline==0) = 100000;
    z_score_all = (mean_fr_bump - mean_fr_baseline)./std_fr_baseline;
    
    
    % get mean angle from z_score data
    mean_fr_ang = zeros(num_units,1);
    max_fr_ang = zeros(num_units,1);
    for i_unit = 1:num_units
        mean_fr_ang(i_unit) = circ_mean(bump_dirs'*pi/180,z_score_all(:,i_unit));
        [~,max_fr_ang(i_unit)] = max(z_score_all(:,i_unit));
    end
    
    % plot move and bump PDs against each other. Plot mod depth for both
    figure()
    subplot(1,3,1)
    plot(pd_move.velPD,pd_bump.velPD,'.')
    subplot(1,3,2)
    plot(pd_move.velModdepth,pd_bump.velModdepth,'.')
    subplot(1,3,3)
    plot(pd_bump.velPD,mean_fr_ang,'.')

    
%% plot z-score for each bump dir for each unit
    for i = 1:size(z_score_all,2)
        figure();
        plot(z_score_all(:,i))
    end
    


%% downsample electrodes to get a uniform distribution of electrodes? 
    % probability of selecting an electrode based on number of electrodes
    % with similar PDs or mean_ang around
    % determine probability by using a circular normal distribution (von
    % mises)
    use_pd = 0;
    percent_elec_keep = 0.6;
    n_groups = 1;
    min_z_score_to_include = 2.0;
    
    if(use_pd)
        data_all = pd_bump.velPD;
    else
        data_all = mean_fr_ang;
    end
    
    include_elec = any(z_score_all > min_z_score_to_include);
    z_score_all_include = z_score_all(:,include_elec);
    data_all = data_all(include_elec);
    chan_list_all_include = chan_list_all(include_elec);
    
    prob_list = zeros(size(data_all));
    num_elec_keep = ceil(numel(prob_list)*percent_elec_keep);
    thetahat=0;
    kappa=10;
    figure(); 
    subplot(1,2,1)
    for i_elec = 1:numel(prob_list)
        data_dist = circ_dist(data_all(i_elec),data_all);
        data_dist_pdf = circ_vmpdf(data_dist,thetahat,kappa);
        prob_list(i_elec) = 1/sum(data_dist_pdf);
        polarplot([0,data_all(i_elec)],[0,prob_list(i_elec)]); hold on
    end

    % go through and select electrodes based on probability
    [data,sample_idx] = datasample(data_all,num_elec_keep,'Weights',prob_list,'Replace',false);

    subplot(1,2,2)
    for i_elec = 1:numel(data)
        polarplot([0,data(i_elec)],[0,1]); hold on
    end
    
    chan_list = chan_list_all_include(sample_idx);
    z_score = z_score_all_include(:,sample_idx);
% sample electrodes based on PDs or mean_fr_ang
    
    
    sample_data = stratifySamplePDs(data,n_groups);
    
%  plot PD around a polar plot, color = repetition idx
    figure(); 
    subplot(1,2,2)
    for i_group = 1:n_groups
        theta_list = data(sample_data.groups(:,i_group));
        for i_plot = 1:numel(theta_list)
            polarplot([0,theta_list(i_plot)],[0,1],'color',getColorFromList(1,i_group-1)); hold on
        end
    end

    % plot PD around a polar plot, color = group_idx
    subplot(1,2,1)
    for i_group_idx = 1:floor(numel(data)/n_groups)
        theta_list = data(sample_data.sort_idx(sample_data.group_idx==i_group_idx));
        for i_plot = 1:numel(theta_list)
            polarplot([0,theta_list(i_plot)],[0,1],'color',getColorFromList(1,i_group_idx-1)); hold on
        end
    end

% plot PD magnitude in all directions for each group
    thetas = linspace(0,2*pi,100);
    neural_mag = zeros(numel(thetas),size(sample_data.groups,2));
    thetas_test = [0,pi/2,pi,3*pi/2];
    best_theta = zeros(size(sample_data.groups,2),1); % overriden
    best_proj = -10000 + zeros(size(best_theta)); % overriden
    figure();
    for i_group = 1:size(neural_mag,2)
        PD_group = data(sample_data.groups(:,i_group));
        proj_data = cos(PD_group-thetas);
        proj_data(proj_data < 0) = 0;
        neural_mag(:,i_group) = mean(proj_data,1,'omitnan')';
        polarplot(thetas,neural_mag(:,i_group),'color',getColorFromList(1,i_group-1),'linewidth',2); hold on
        % get best set of 4 directions separated by pi/2
        for i_theta = 1:numel(thetas)
            thetas_ = thetas_test + thetas(i_theta);
            proj_data = cos(PD_group-thetas_);
            proj_data(proj_data < 0) = nan;
            theta_proj = mean(mean(proj_data,'omitnan'));
            if(theta_proj > best_proj(i_group))
                best_proj(i_group) = theta_proj;
                best_theta(i_group) = thetas(i_theta);
            end
        end
%         for i_plot = 1:numel(thetas_test)
%             polarplot([0,best_theta(i_group)+thetas_test(i_plot)],[0,1],'linestyle','--','linewidth',1,'color',getColorFromList(1,i_group-1));
%         end
    end

    
%% build patterns
    constant_scramble = 0;
    tgt_dirs = (45+[0,90,180,270])*pi/180;
    group_idx = 1;
    min_z_score = 0;
    max_z_score = 1;
        
    min_stim_norm = 0.1; max_stim_norm = 1;
    
    num_scrambles = 100000;
    min_MSD_percentile = 75; % percent
    
    % set up pattern data struct
    pattern_data = [];
    pattern_data.num_patterns = numel(tgt_dirs)*2; % must be even (pairs of bio and nonbio)
        
    pattern_data.pred_dirs = [];
    pattern_data.is_biomimetic = [];
    pattern_data.pattern = [];
    
    % setup heatmap inputs
    heatmap_input_data = [];
    heatmap_input_data.map_filename = input_data.mapFileName;
    heatmap_input_data.num_colors = 100;
    heatmap_input_data.base_color = [153,51,255]/256; % max stim value is this color
    
    % normalize z_score for each neuron
    z_score_norm = 1.5*(z_score - min(z_score))./(max(z_score) - min(z_score)) - 0.5;
    
    % make pattern_data, biomimetic patterns first
    if(use_pd)
        % for each tgt dir, project PD onto tgt_dir, if < 0, set to 0
        for i_tgt = 1:numel(tgt_dirs)
            PD_proj = cos(tgt_dirs(i_tgt) - data(sample_data.groups(:,group_idx)));
            PD_proj(PD_proj < 0) = 0;
            
            pattern_data.pred_dirs(end+1) = tgt_dirs(i_tgt);
            bio_pattern = []; 
            bio_pattern.stim_norm = PD_proj; 
            bio_pattern.chans = chan_list(sample_data.groups(:,group_idx),1); % convert idx into channel number
            pattern_data.pattern{end+1} = bio_pattern;
            pattern_data.is_biomimetic(end+1) = 1;
        end
    else % use delta FR during bumps
        for i_tgt = 1:numel(tgt_dirs)
            % get idx in z_score
            [~,z_score_idx] = min(abs(circ_dist(tgt_dirs(i_tgt),bump_dirs*pi/180)));
            z_score_proj = z_score_norm(z_score_idx,sample_data.groups(:,group_idx))';
            z_score_proj(z_score_proj < min_z_score) = 0;
            % set 0-1 for stim_norm
            stim_norm = z_score_proj; %(z_score_proj - min_z_score)/(max_z_score - min_z_score);
            stim_norm(stim_norm < 0) = 0;
            stim_norm(stim_norm < min_stim_norm & stim_norm > 0) = min_stim_norm;
            stim_norm(stim_norm > max_stim_norm) = max_stim_norm;
            
            if(sum(stim_norm >0) > 16)
                % remove electrodes
                num_remove = sum(stim_norm>0) - 16;
                non_zero_idx = find(stim_norm>0);
                [~,make_zero] = datasample(non_zero_idx,num_remove,'Replace',false);
                stim_norm(non_zero_idx(make_zero)) = 0;
                disp(num_remove);
            end
            
            pattern_data.pred_dirs(end+1) = tgt_dirs(i_tgt);
            bio_pattern = []; 
            bio_pattern.stim_norm = stim_norm; 
            bio_pattern.chans = chan_list(sample_data.groups(:,group_idx),1);
           
            pattern_data.pattern{end+1} = bio_pattern;
            pattern_data.is_biomimetic(end+1) = 1;
        end
    end
    % make nonbiomimetic patterns by scrambling biomimetic patterns
    % first scramble biomimetic patterns a lot to generate a distribution
    % of mean squared differences
    min_MSD_list = zeros(num_scrambles,1);
        
    for i_scramble = 1:num_scrambles
        pattern_idx = ceil(rand()*numel(pattern_data.pattern));
        scramble_idx = randperm(numel(pattern_data.pattern{pattern_idx}.chans));
        
        scramble_stim_norm = pattern_data.pattern{pattern_idx}.stim_norm(scramble_idx);
        
        min_MSD = 1000000;
        for i_pattern = 1:numel(pattern_data.pattern)
            MSD = mean((scramble_stim_norm-pattern_data.pattern{i_pattern}.stim_norm).^2);
            if(MSD < min_MSD)
                min_MSD = MSD;
            end
        end
        min_MSD_list(i_scramble) = min_MSD;
    end
    threshold_min_MSD = prctile(min_MSD_list,min_MSD_percentile);
    
    
    % then shuffle each biomimetic pattern until each min MSD is more than
    % the threshold
    
    if(constant_scramble) % scramble all patterns
        min_MSD = threshold_min_MSD -1; % to get in to the loop
        pattern_counter = 0;
        best_MSD = 0;
        best_scramble_idx = 0;
        while(min_MSD < threshold_min_MSD && pattern_counter < 100000)
            scramble_idx = randperm(size(sample_data.groups,1));
            min_MSD = 10000000; 
            
            for i_pattern = 1:sum(pattern_data.is_biomimetic)
                scramble_stim_norm = pattern_data.pattern{pattern_idx}.stim_norm(scramble_idx);
                for i_temp = 1:sum(pattern_data.is_biomimetic)
                    MSD = mean((scramble_stim_norm-pattern_data.pattern{i_temp}.stim_norm).^2);
                    if(MSD < min_MSD)
                        min_MSD = MSD;
                    end
                end
            end
            if(min_MSD >= best_MSD)
                best_scramble_idx = scramble_idx;
                best_MSD = min_MSD;
            end
            pattern_counter = pattern_counter + 1;
        end
        disp(pattern_counter)
        for i_pattern = 1:sum(pattern_data.is_biomimetic)
            pattern_data.pattern{end+1}.stim_norm = pattern_data.pattern{i_pattern}.stim_norm(best_scramble_idx);
            pattern_data.pattern{end}.chans = pattern_data.pattern{i_pattern}.chans;
            pattern_data.is_biomimetic(end+1) = 0;
            pattern_data.pred_dirs(end+1) = pattern_data.pred_dirs(i_pattern);
        end
    else % scramble each pattern individually
        for i_pattern = 1:sum(pattern_data.is_biomimetic)
            min_MSD = threshold_min_MSD -1; % to get in to the loop
            pattern_counter = 0;
            best_MSD = 0;
            best_scramble_idx = 0;
            while(min_MSD < threshold_min_MSD && pattern_counter < 100000)
                min_MSD = 10000000;
                scramble_idx = randperm(numel(pattern_data.pattern{pattern_idx}.chans));
                scramble_stim_norm = pattern_data.pattern{pattern_idx}.stim_norm(scramble_idx);
                for i_temp = 1:sum(pattern_data.is_biomimetic)
                    MSD = mean((scramble_stim_norm-pattern_data.pattern{i_temp}.stim_norm).^2);
                    if(MSD < min_MSD)
                        min_MSD = MSD;
                    end
                end
                if(min_MSD >= best_MSD)
                    best_scramble_idx = scramble_idx;
                    best_MSD = min_MSD;
                end
                pattern_counter = pattern_counter + 1;
            end
            pattern_data.pattern{end+1}.stim_norm = pattern_data.pattern{i_pattern}.stim_norm(best_scramble_idx);
            pattern_data.pattern{end}.chans = pattern_data.pattern{i_pattern}.chans;
            pattern_data.is_biomimetic(end+1) = 0;
            pattern_data.pred_dirs(end+1) = pattern_data.pred_dirs(i_pattern);
            disp(pattern_counter)
        end
    end

%% plot patterns
    plotPatternDataHeatmap(pattern_data, heatmap_input_data);

    % get pattern metrics
    num_patterns = numel(pattern_data.pattern);
    total_stim = zeros(num_patterns,1);
    num_stim_chans = zeros(num_patterns,1);
    predicted_dir_mag = zeros(num_patterns,3);
    predicted_dir_mag(:,1) = pattern_data.pred_dirs;
    for i_pattern = 1:num_patterns
        % total "stimulation" defined as sum of stim_norm
        total_stim(i_pattern) = sum(pattern_data.pattern{i_pattern}.stim_norm);
        % # channels actually stimulated
        num_stim_chans(i_pattern) = sum(pattern_data.pattern{i_pattern}.stim_norm > 0);
        % predicted dir -- mean of angle weighted by stim_norm
        predicted_dir_mag(i_pattern,2) = circ_mean(data(sample_data.groups(:,group_idx)),pattern_data.pattern{i_pattern}.stim_norm);
        % predicted mag -- sum of stim_norm projected onto
            % predicted_tgt_dir
        predicted_dir_mag(i_pattern,3) = sum(pattern_data.pattern{i_pattern}.stim_norm.*cos(data(sample_data.groups(:,group_idx))-predicted_dir_mag(i_pattern,1)));
    end
    
    disp(num_stim_chans)
    disp(total_stim)
    disp(predicted_dir_mag)
    
%% convert pattern data to stim_array for the cerestim
    stim_params = [];
    stim_params.IPI = 3/1000; % s
    stim_params.max_freq = 300; % Hz
    stim_params.min_freq = 0; % Hz
    stim_params.train_length = 0.5; % s
    stim_array_data = makeStimArrayWrapper(pattern_data,stim_params);

    % make stim_array which is formatted properly
    stim_array = cell(numel(stim_array_data.stim_array),1);
    for i_pattern = 1:numel(stim_array)
        stim_array{i_pattern}.stim_pattern = stim_array_data.stim_array{i_pattern};
        stim_array{i_pattern}.chans = stim_array_data.chans{i_pattern};
    end

%% plot mean IPI and list of IPIs with desired IPI for each electrode
    figure();
    for i_pattern = 1:numel(stim_array_data.stim_array)
        subplot(4,2,i_pattern); hold on
        for i_elec = 1:size(stim_array_data.stim_array{i_pattern},1)
            IPI_list = diff(find(stim_array_data.stim_array{i_pattern}(i_elec,:)))*stim_params.IPI;
            if(~isempty(IPI_list))
                mean_IPI = mean(IPI_list);
                IPI_list = unique(IPI_list);

                plot(i_elec,mean_IPI,'r.','markersize',12);
                plot(i_elec,IPI_list,'bo','markersize',8);
            end
            
            % get desired IPI
            chan_idx = find(pattern_data.pattern{i_pattern}.chans == stim_array_data.chans{i_pattern}(i_elec));
            desired_IPI = 1/(pattern_data.pattern{i_pattern}.stim_norm(chan_idx)*stim_params.max_freq);
            plot(i_elec,desired_IPI,'k*','markersize',8);
        end
    end

    
%% convert pattern data to wave_mappings
% wave_mapping{code}(1) = channel, wave_mapping{code}(2) = desired normalized stim frew,
% wave_mapping{code}(3) = freq code for cerestim (1-15)

    num_freqs = 15;
    freq_all = linspace(0,1,num_freqs+1); freq_all = freq_all(2:end);
    wave_mappings = {};
    
    
    for i_patt = 1:numel(pattern_data.pattern)
        keep_mask = pattern_data.pattern{i_patt}.stim_norm > 0;
        stim_norm = pattern_data.pattern{i_patt}.stim_norm(keep_mask);
        chans = pattern_data.pattern{i_patt}.chans(keep_mask);
        
        % round stim_norm to values in freq_all
        stim_norm_adj = zeros(numel(stim_norm),1);
        freq_code = zeros(numel(stim_norm),1);
        
        for i_chan = 1:numel(stim_norm_adj)
            [~,min_idx] = min(abs(stim_norm(i_chan) - freq_all)); 
            stim_norm_adj(i_chan) = freq_all(min_idx);
            freq_code(i_chan) = min_idx;
        end
        
        wave_mappings{i_patt} = [chans,stim_norm,freq_code];
        
    end

    
    