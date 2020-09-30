%% set filename
    input_data.folderpath = 'D:\Lab\Data\LFP\Han_20191101_CObumpmove\';
    
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    
    input_data.date = '20191106';
    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyDuncan';
    input_data.ranBy = 'ranByJoe';
    input_data.lab = 6;  
    input_data.mapFile = strcat('mapFile',mapFileName);
    input_data.task = 'taskCObump';

    pwd=cd;
    input_data.fileNum = 1;
  
    input_data.center_y = -33;
    input_data.center_x = 3;
    input_data.num_bins = 8;
    
%% make cds and trial data
   % REMOVE ID FROM UNITS and make trial data

    cd(input_data.folderpath)

    file_name = dir('*nev*');

    params.event_list = {'startTime';'goCueTime';'bumpTime';'bumpDir'};
    params.trial_results = {'R'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [255];

    move_onset_params.pre_move_thresh = 1000;
    move_onset_params.min_s = 3;
    move_onset_params.start_idx_offset = -10;
    move_onset_params.max_rt_offset = 400;
    
    cds = commonDataStructure();
    cds.file2cds(strcat(input_data.folderpath,file_name(1).name),input_data.array,input_data.monkey,input_data.ranBy,...
        input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
   
    td_all = parseFileByTrial(cds,params);
    td_all = stripSpikeSorting(td_all);
    td_all = getSpeed(td_all);
    td_all = removeBadTrials(td_all);
    
    if(td_all(1).bin_size < 0.05)
        % set it to 50ms
        td_all = binTD(td_all,ceil(0.05/td_all(1).bin_size));
    end
  
    
    td_all = getMoveOnset(td_all,move_onset_params);
    td_all = removeBadTrials(td_all);


%% preprocess lfp to go with trial data

    % subtract regressed common average from each electrode
    % FFT in 256ms (?) windows
    % calculate power in frequency bands 
    input_data.subtract_common_average = 1;
    input_data.FFT_length = 256; % in ms
    input_data.power_bands = [15,30];
    
    [td_lfp, lfp_data] = preProcessLFP(cds, td_all ,input_data);
    
 
%% plot power data and neural firing during reaches
    plot_data = 1;
    power_colors = inferno(size(input_data.power_bands,1)+1);
    power_idx = 1;
    reach_dir = unique([td_lfp.target_direction]);
    window = [-0.5,1.5]; % around go cue
    spikes_subplot_idx = [15,3,11,23];
    power_subplot_idx = [14,8,12,18];
    
    mean_firing_rate = zeros(size(td_lfp(1).LeftS1_unit_guide,1),numel(reach_dir),floor(diff(window)/td_lfp(1).bin_size));
    mean_power = zeros(size(td_lfp(1).LeftS1_unit_guide,1),numel(reach_dir),size(input_data.power_bands,1), floor(diff(window)/td_lfp(1).bin_size));
    x_data = linspace(window(1),window(2),floor(diff(window)/td_lfp(1).bin_size));
    for i_unit = 21%:size(td_all(1).LeftS1_unit_guide,1)
        if(plot_data) figure('Position',[680 266 856 712]); end
        for i_reach = 1:numel(reach_dir)
            trial_idx = find([td_lfp.target_direction] == reach_dir(i_reach) & [td_lfp.result] == 'R' & isnan([td_lfp.bumpDir]));
            for trial = trial_idx
                trial_window = td_lfp(trial).idx_goCueTime + floor(window/td_lfp(1).bin_size);
                trial_window(end) = trial_window(end)-1;
                % get mean firing rate (add up spike counts, divide after
                % for loop)
                mean_firing_rate(i_unit,i_reach,:) = mean_firing_rate(i_unit,i_reach,:) + ...
                    reshape(td_lfp(trial).LeftS1_spikes(trial_window(1):trial_window(2),i_unit),1,1,size(mean_firing_rate,3));
            
                % get mean power (add up power, divide after for loop)
                mean_power(i_unit,i_reach,:,:) = squeeze(mean_power(i_unit,i_reach,:,:)) + ...
                    squeeze(td_lfp(trial).lfp_data(trial_window(1):trial_window(2),i_unit,:));
            end
            
            % normalize by number of trials and bin size
            mean_firing_rate(i_unit,i_reach,:) = mean_firing_rate(i_unit,i_reach,:)/numel(trial_idx)/td_lfp(1).bin_size;
            % normalzie power by number of trials
            mean_power(i_unit,i_reach,:,:) = mean_power(i_unit,i_reach,:,:)/numel(trial_idx);
            
            
            if(plot_data)
                % plot firing rate
                ax = subplot(5,5,spikes_subplot_idx(i_reach));
                plot(x_data,squeeze(mean_firing_rate(i_unit,i_reach,:)),'color',getColorFromList(1,1),'linewidth',1.5);
                hold on;
                plot([0,0],ax.YLim,'r--','linewidth',0.75)
                xlim(window)
                % plot power
                ax = subplot(5,5,power_subplot_idx(i_reach));
                for i_power = power_idx%:size(mean_power,3)
                    plot(x_data,squeeze(mean_power(i_unit,i_reach,i_power,:)),'color',power_colors(i_power,:),'linewidth',1.5); hold on;
                end
                plot([0,0],ax.YLim,'r--','linewidth',0.75)
                xlim(window)
            end
        end
    end

    
%% try to decode hand kinematics (velocity and position) with spikes and LFP to see if LFP data is even good....
    percent_train = 0.8;
    num_lags = 10;
    window = [-0.5,1.5]; % relative to go cue
    window_idx = floor(window/td_lfp(1).bin_size);
    
    train_idx = datasample((1:1:numel(td_lfp)),floor(numel(td_lfp)*percent_train),'Replace',false);
    
    % concatenate trials (training and testing)
    
    train_mask = [];
    pos = [];
    vel = [];
    lfp = [];
    spike_count = [];
    
    for i_trial = 1:numel(td_lfp)
        trial_window = td_lfp(i_trial).idx_goCueTime + window_idx;
        temp_spike = [];
        temp_lfp = [];
        for i_lag = 0:num_lags-1
            temp_spike = [temp_spike, td_lfp(i_trial).LeftS1_spikes(trial_window(1)-i_lag:trial_window(2)-i_lag,:)];
            temp_lfp = [temp_lfp, reshape(td_lfp(i_trial).lfp_data(trial_window(1)-i_lag:trial_window(2)-i_lag,:,:),...
                diff(trial_window)+1,size(td_lfp(i_trial).lfp_data,2)*size(td_lfp(i_trial).lfp_data,3))];
        end
        
        pos = [pos;td_lfp(i_trial).pos(trial_window(1)-i_lag:trial_window(2)-i_lag,:)];
        vel = [vel;td_lfp(i_trial).vel(trial_window(1)-i_lag:trial_window(2)-i_lag,:)];
        spike_count = [spike_count;temp_spike];
        lfp = [lfp; temp_lfp];
        train_mask = [train_mask; any(train_idx==i_trial)+zeros(diff(trial_window)+1,1)];
    end

    [lfp_pca_coeff,lfp,~,~,explained] = pca(lfp,'NumComponents',960);
    [spike_pca_coeff,spike_count] = pca(spike_count,'NumComponents',960);
    
%% build glm

    y_data = [pos,vel]; test_pred = y_data(train_mask==0,:);
    lfp_mdl = {}; spike_mdl = {};
    lfp_pred = []; spike_pred = [];
    lfp_test_r2 = []; spike_test_r2 = [];
    for i_var = 1:4
        lfp_mdl{i_var} = fitrlinear(lfp(train_mask==1,:),y_data(train_mask==1,i_var),'Learner','leastsquares');
        lfp_pred(:,i_var) = predict(lfp_mdl{i_var},lfp(train_mask==0,:));
        lfp_test_r2(i_var) = 1-(sum((lfp_pred(:,i_var)-test_pred(:,i_var)).^2))/(sum((test_pred(:,i_var)-mean(test_pred(:,i_var))).^2));
        
        spike_mdl{i_var} = fitrlinear(spike_count(train_mask==1,:),y_data(train_mask==1,i_var),'Learner','leastsquares');
        spike_pred(:,i_var) = predict(spike_mdl{i_var},spike_count(train_mask==0,:));
        spike_test_r2(i_var) = 1-(sum((spike_pred(:,i_var)-(test_pred(:,i_var))).^2))/(sum((test_pred(:,i_var)-mean(test_pred(:,i_var))).^2));
    end
    
    % test glm
    
    

%% consistency of LFP and neural discharge during trials (cds)
    trial_idx = find(isnan([cds.trials.bumpDir]) & [cds.trials.tgtDir] == 90 & [cds.trials.result == 'R']);
    trial_idx = trial_idx(1:10);
    elec_idx = 33;
    
    figure();
    for i_trial = 1:numel(trial_idx)
        subplot(numel(trial_idx),1,i_trial)
        window = cds.trials.goCueTime(trial_idx(i_trial)) + [-1,2];

        lfp_idx = [find(cds.lfp.t > window(1),1,'first'),find(cds.lfp.t > window(2),1,'first')];

        spike_mask = cds.units(elec_idx).spikes.ts > window(1) & cds.units(elec_idx).spikes.ts < window(2);
        spike_ts = cds.units(elec_idx).spikes.ts(spike_mask) - cds.trials.goCueTime(trial_idx(i_trial));

        plot(cds.lfp.t(lfp_idx(1):lfp_idx(2))-cds.trials.goCueTime(trial_idx(i_trial)),lfp_data(lfp_idx(1):lfp_idx(2),elec_idx))
        hold on

        plot([spike_ts,spike_ts]',(max(lfp_data(lfp_idx(1):lfp_idx(2),elec_idx))+repmat([1,11],numel(spike_ts),1))','k')

        xlim(window-cds.trials.goCueTime(trial_idx(i_trial)))
    end

%% get LFP on all trials around go cue
    trial_idx = find(isnan([cds.trials.bumpDir]) & [cds.trials.tgtDir] == 0 & [cds.trials.result == 'R']);
    lfp_go_cue_data = [];
    for i_trial = trial_idx'
        window = cds.trials.goCueTime(i_trial) + [-1,2];
        lfp_idx = [find(cds.lfp.t > window(1),1,'first'),find(cds.lfp.t > window(2),1,'first')];
        lfp_go_cue_data(end+1,:,:) = lfp_data(lfp_idx(1):lfp_idx(2),:);
    end
    go_cue_idx = find(cds.lfp.t > cds.trials.goCueTime(i_trial),1,'first') - lfp_idx(1);

    
    
    
    