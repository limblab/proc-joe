%% set initial parameters

    input_data.folderpath = 'D:\Lab\Data\StimPDs\Han_20191003_CObump_stimDuringTask\';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\6250-002085\SN 6250-002085.cmp';
    
    input_data.date = '20191003';
    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyHan';
    input_data.ranBy = 'ranByJoe';
    input_data.lab = 6;
    input_data.mapFile = strcat('mapFile',mapFileName);
    input_data.task = 'taskCObump';

    pwd=cd;
    input_data.fileNum = 1;
  
    input_data.center_y = -33;
    input_data.center_x = 3;
    input_data.num_bins = 8;
    
%% make cds
    cd(input_data.folderpath)

    file_name = dir('*nev*');
    
    params.event_list = {'goCueTime';'bumpTime';'bumpDir';'stimTime';'stimCode'};
    params.trial_results = {'R'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [0,255];
    
    td_all = [];
    td_all = [];
    for file_num = 1:numel(file_name)
        cds = commonDataStructure();
        cds.file2cds(strcat(input_data.folderpath,file_name(file_num).name),input_data.array,input_data.monkey,input_data.ranBy,...
            input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');


        td_temp = parseFileByTrial(cds,params);
        td_temp = getSpeed(td_temp);
        td_all = [td_all,td_temp];

    end 

    td_all = getAdjustedSpikeCount(td_all,struct('num_pulses',10,'IPI',10.6));
    td_all = getNorm(td_all,struct('signals',{'vel'}));
    td_all = getMoveOnsetAndPeak(td_all,struct('which_field','vel_norm','start_idx_offset',-5));
    td_all = removeBadTrials(td_all);
    
    % remove stim trials where stim happens too far after movement onset
    td_all = td_all(isnan([td_all.idx_stimTime]) | ... 
        [td_all.idx_stimTime] - [td_all.idx_movement_on] <= 0.4/td_all(1).bin_size);
    
%% plot kinematics....
    % plot example traces for each target
    % plot velocity for stim and non-stim trials for each target
    tgt_dir = unique([td_all.target_direction]);
    subplot_idx = [6,2,4,8];
    num_trials = 20;
    color_list = getColorFromList(1,0:3);
    
    f_trace = figure(); hold on;
    formatForLee(gcf); xlabel('X-pos (cm)'); ylabel('Y-pos (cm)');
    
    f_vel = figure(); hold on;
    formatForLee(gcf); xlabel('Time relative to movement onset (s)');
    ylabel('Speed (cm/s)');
    for i_tgt = 1:numel(tgt_dir)
        td_dir = td_all([td_all.target_direction] == tgt_dir(i_tgt));
        td_dir = td_dir(datasample(1:1:numel(td_dir),num_trials));
        td_dir = trimTD(td_dir,{'idx_movement_on',-4},{'idx_endTime',0});
        
        for i_trial = 1:numel(td_dir)
%             % plot traces
            figure(f_trace);
            linestyle = '-';
            if(~isnan(td_dir(i_trial).idx_stimTime) && td_dir(i_trial).idx_stimTime > td_dir(i_trial).idx_goCueTime)
                stim_idx = td_dir(i_trial).idx_stimTime;
                plot(td_dir(i_trial).pos(stim_idx,1),td_dir(i_trial).pos(stim_idx,2),'k.','markersize',12)
                linestyle = '--';
            end
            plot(td_dir(i_trial).pos(:,1),td_dir(i_trial).pos(:,2),'color',color_list(i_tgt,:),'linestyle',linestyle);
            
            % plot velocity on correct subplot
            figure(f_vel);
            subplot(3,3,subplot_idx(i_tgt)); hold on;
            if(~isnan(td_dir(i_trial).idx_stimTime) && td_dir(i_trial).idx_stimTime > td_dir(i_trial).idx_goCueTime)
                linestyle = '--';
            else
                linestyle = '-';
            end
            x_data = ((1:1:numel(td_dir(i_trial).speed)) - td_dir(i_trial).idx_movement_on - 1)*td_dir(i_trial).bin_size;
            plot(x_data,td_dir(i_trial).speed,'color',color_list(i_tgt,:),'linestyle',linestyle);
            if(~isnan(td_dir(i_trial).idx_stimTime) && td_dir(i_trial).idx_stimTime > td_dir(i_trial).idx_goCueTime)
                stim_idx = td_dir(i_trial).idx_stimTime;
                plot(x_data(stim_idx),td_dir(i_trial).speed(stim_idx),'k.','markersize',12)
            end
            xlim([-0.2,0.6])
            ylim([0,75]);
            xlabel('time after movement onset (s)');
            ylabel('speed (cm/s)');
        end
        
    end
    
        
%% get PSTH data for each target dir and stim trials
    num_pulses = 10;
    stim_length = ceil(((num_pulses-1)*10.6/1000)/td_all(1).bin_size);
    
    tgt_dir = unique([td_all.target_direction]);
    array_name = input_data.array(6:end);
    num_neurons = size(td_all(1).([array_name,'_unit_guide']),1);
    
    data_window = ceil([-0.25,0.8]/td_all(1).bin_size); % window about alignment to get FR across trials
    
    mean_fr_move = zeros(diff(data_window)+1,num_neurons,numel(tgt_dir),2); % time x neurons x tgts x no stim/stim
    mean_speed_move = zeros(diff(data_window)+1,numel(tgt_dir),2); % time x tgts x no stim/stim
    mean_fr_hold = zeros(diff(data_window)+1,num_neurons,2); % time x neurons x no stim/stim
    num_hold_trials = sum(~isnan([td_all.idx_stimTime]) | [td_all.idx_stimTime] < [td_all.idx_goCueTime]);
    num_move_trials = sum(~isnan([td_all.idx_stimTime]) | [td_all.idx_stimTime] > [td_all.idx_goCueTime]);
    fr_hold_all = zeros(num_neurons,num_hold_trials); % mean across trials
    mean_speed_hold = zeros(diff(data_window)+1,2); % time  x no stim/stim
    
    move_fr_all = zeros(num_neurons,num_move_trials);
    move_fr_tgt_dir = zeros(num_move_trials,1);
    move_fr_stim_mask = zeros(num_move_trials,1);
    
    chan_num = td_all(1).([array_name,'_unit_guide'])(:,1);
    num_trials_move = zeros(numel(tgt_dir),2); num_trials_hold = [0,0];
    
    % for non-stim hold trials, plot time relative to when stim typically occurs
    stim_hold_idx = [td_all.idx_goCueTime] - [td_all.idx_stimTime];
    move_to_stim_idx = abs(floor(mean(stim_hold_idx(~isnan(stim_hold_idx)& stim_hold_idx < 0))));
    hold_idx = floor(mean(stim_hold_idx(~isnan(stim_hold_idx) & stim_hold_idx > 0))); % window relative to this before go cue for non stim trials
    fr_trial_counter = 1;
    move_trial_counter = 1;
    
    hold_stim_mask = zeros(num_hold_trials,1);
    move_stim_mask = zeros(num_move_trials,1);
    
    is_hold_stim = [];
    
    for i_trial = 1:numel(td_all)
        stim_idx = ~isnan(td_all(i_trial).idx_stimTime) + 1; % 1 = no stim, 2 = stim
        is_hold_stim = stim_idx == 2 & td_all(i_trial).idx_movement_on > td_all(i_trial).idx_stimTime;
        
        % deal with hold data. if a stim trial, window relative to stim
        % timing. Otherwise use hold_idx relative to go_cue
        switch stim_idx
            case 1 % no stim
                idx_middle = td_all(i_trial).idx_goCueTime - hold_idx - move_to_stim_idx;
            case 2 % stim
                idx_middle = td_all(i_trial).idx_stimTime - move_to_stim_idx + 1;
        end
        if(stim_idx == 1 || is_hold_stim)
            num_trials_hold(stim_idx) = num_trials_hold(stim_idx)+1;
            mean_fr_hold(:,:,stim_idx) = mean_fr_hold(:,:,stim_idx) + ...
                td_all(i_trial).([array_name,'_spikes'])(idx_middle+data_window(1):idx_middle+data_window(2),:)/td_all(i_trial).bin_size;
            mean_speed_hold(:,stim_idx) = mean_speed_hold(:,stim_idx) + ...
                td_all(i_trial).speed(idx_middle+data_window(1):idx_middle+data_window(2),:);
            
            fr_hold_all(:,fr_trial_counter) = mean(td_all(i_trial).([array_name,'_spikes'])(idx_middle+move_to_stim_idx-1:idx_middle+move_to_stim_idx-1+stim_length-1,:)/td_all(i_trial).bin_size);
            
            hold_stim_mask(fr_trial_counter) = stim_idx == 2;
            fr_trial_counter = fr_trial_counter + 1;
        end
%         % deal with move data
%         % for move data, get time relatve to go cue for both stim and non
%         % stim trials
        idx_middle = td_all(i_trial).idx_goCueTime;
        tgt_idx = find(td_all(i_trial).target_direction == tgt_dir);
        num_trials_move(tgt_idx,stim_idx) = num_trials_move(tgt_idx,stim_idx) + 1;
        mean_fr_move(:,:,tgt_idx,stim_idx) = mean_fr_move(:,:,tgt_idx,stim_idx) + ...
            td_all(i_trial).([array_name,'_spikes'])(idx_middle+data_window(1):idx_middle+data_window(2),:)/td_all(i_trial).bin_size;
        mean_speed_move(:,tgt_idx,stim_idx) = mean_speed_move(:,tgt_idx,stim_idx) + ...
            td_all(i_trial).speed(idx_middle+data_window(1):idx_middle+data_window(2),:);
        if(stim_idx == 2) % use stim time
            idx_middle = td_all(i_trial).idx_stimTime;
            move_fr_all(:,move_trial_counter) = mean(td_all(i_trial).([array_name,'_spikes'])(idx_middle:idx_middle+stim_length-1,:)/td_all(i_trial).bin_size);
        else % add mean stim time to go cue
            move_fr_all(:,move_trial_counter) = mean(td_all(i_trial).([array_name,'_spikes'])(idx_middle+move_to_stim_idx:idx_middle+move_to_stim_idx+stim_length-1,:)/td_all(i_trial).bin_size);
        end
            move_fr_tgt_dir(move_trial_counter) = tgt_idx;
        move_fr_stim_mask(move_trial_counter) = stim_idx == 2;
        move_trial_counter = move_trial_counter + 1;
    end
    
    % normalize by num_trials
    
    mean_fr_move = mean_fr_move./repmat(reshape(num_trials_move,1,1,size(num_trials_move,1),size(num_trials_move,2)),size(mean_fr_move,1),size(mean_fr_move,2));
    mean_fr_hold = mean_fr_hold./repmat(reshape(num_trials_hold,1,size(num_trials_hold,1),size(num_trials_hold,2)),size(mean_fr_hold,1),size(mean_fr_hold,2));
    
    mean_speed_move = mean_speed_move./repmat(reshape(num_trials_move,1,size(num_trials_move,1),size(num_trials_move,2)),size(mean_speed_move,1),1,1);
    mean_speed_hold = mean_speed_hold./repmat(reshape(num_trials_hold,1,numel(num_trials_hold)),size(mean_speed_hold,1),1);

    
%% plot raster for each neuron
    num_neurons = size(td_all(1).LeftS1_spikes,2);
    
    spike_row = []; spike_time = [];
    stim_row = []; stim_time = [];
    
    for i_unit = 1:num_neurons


    end
%% plot PSTH data for each neuron
    move_subplot_idx = [6,2,4,8]; % right, up, left, down
    hold_subplot_idx = 5;
    
    for i_unit = 1:num_neurons
        ax_list = {};
        f=figure('Position',[2164 213 905 695]);
        f.Name = [input_data.monkey(7:end),'_cobump_chan25stim_chan',num2str(chan_num(i_unit)),'rec_PSTH'];
        % plot 5 figures, one for each tgt dir and one for center hold
        % plot stim and no stim on each figure
        
        % plot hold data and stim time
        x_data = (round(data_window(1)):round(data_window(end)))*td_all(1).bin_size;
        x_stim = move_to_stim_idx*td_all(1).bin_size - td_all(1).bin_size/2;
        
        ax_list{1}=subplot(3,3,hold_subplot_idx); hold on
        plot(x_data,mean_fr_hold(:,i_unit,1),'k');
        plot(x_data,mean_fr_hold(:,i_unit,2),'r');
        
        y_lim_use = ax_list{1}.YLim;
        plot([x_stim,x_stim],[-1000,1000],'r--')
        
        % clean up plot
        ylabel('Firing rate (Hz)')
        xlabel('Time during center hold (s)');
        formatForLee(gcf);
        
        for i_tgt = 1:numel(tgt_dir)
            ax_list{i_tgt+1}=subplot(3,3,move_subplot_idx(i_tgt)); hold on
            plot(x_data,mean_fr_move(:,i_unit,i_tgt,1),'k');
            plot(x_data,mean_fr_move(:,i_unit,i_tgt,2),'r');

            xlim([x_data(1),x_data(end)])
            
            % update y_lim_use
            y_lim_use = [min(y_lim_use(1),ax_list{i_tgt+1}.YLim(1)),max(y_lim_use(2),ax_list{i_tgt+1}.YLim(2))];
            
            % plot stim data
            plot([x_stim,x_stim],[-1000,1000],'r--')
            
            % clean up plot
            xlabel('Time after movement onset (s)')
            ylabel('Firing rate (Hz)')
            formatForLee(gca);
        end
        
        linkaxes([ax_list{:}],'xy');
        ylim(y_lim_use);
        xlim([x_data(1),x_data(end)])
        
    end

%% plot firing rate during stim for center hold vs movement
    
    stim_bins = abs(data_window(1)) + 1 + (move_to_stim_idx:1:(move_to_stim_idx+stim_length));
    stim_chan = 25;
    alpha = 0.05;
    subtract_move_fr = 0;
    
    center_hold_fr = []; center_hold_std = []; move_fr_center_fr_diff = [];
    move_fr = []; move_std = []; mod_depth = []; move_fr_baseline = [];
    is_responsive = [];
    stim_chan_mask = [];
    for i_unit = 1:num_neurons
        if(ranksum(fr_hold_all(i_unit,hold_stim_mask==1),fr_hold_all(i_unit,hold_stim_mask==0),'tail','right') < 0.01)
            is_responsive(i_unit) = 1;
            sub_val = 0;
            if(subtract_move_fr)
                sub_val = mean(fr_hold_all(i_unit,hold_stim_mask==0));
            end
            
            center_hold_fr(end+1:end+numel(tgt_dir),1) = mean(fr_hold_all(i_unit,hold_stim_mask==1) - sub_val);
            center_hold_std(end+1:end+numel(tgt_dir),1) = std(fr_hold_all(i_unit,hold_stim_mask==1) - sub_val);

            
            % max dir
            for i_tgt = 1:numel(tgt_dir)
                if(subtract_move_fr)
                    sub_val = mean(move_fr_all(i_unit,move_fr_stim_mask==0 & move_fr_tgt_dir==i_tgt));
                end
                move_fr(end+1,1) = mean(move_fr_all(i_unit,move_fr_stim_mask==1 & move_fr_tgt_dir==i_tgt)) - sub_val;
                move_fr_baseline(end+1,1) = mean(move_fr_all(i_unit,move_fr_stim_mask==0 & move_fr_tgt_dir==i_tgt));
                move_std(end+1,1) = std(move_fr_all(i_unit,move_fr_stim_mask==1 & move_fr_tgt_dir==i_tgt));
                mod_depth(end+1,1) = mean(move_fr_all(i_unit,move_fr_stim_mask==0 & move_fr_tgt_dir==i_tgt)) - mean(fr_hold_all(i_unit,hold_stim_mask==0));
            
                stim_chan_mask(end+1,1) = chan_num(i_unit) == stim_chan;
                
                move_fr_center_fr_diff(end+1,1) = move_fr(end) - center_hold_fr(end);
            end
        end
    end
    
    
    figure(); hold on
    ax = gca;
    num_colors = 200;
    inferno_list  = inferno(num_colors+20);
    color_idx = floor((num_colors-1)*(mod_depth - min(mod_depth))/(max(mod_depth)-min(mod_depth))) + 1;
    color_list = inferno_list(color_idx,:);
    ax.ColorOrder = color_list;
    ax.ColorOrderIndex = 1;
    
    for i = 1:numel(move_fr)
        if(stim_chan_mask(i)==1)
            marker = '*'; markersize = 8;
        else
            marker = '.'; markersize = 12;
        end
        plot(center_hold_fr(i),move_fr(i),...
            marker,'markersize',markersize);
    end

    tick_labels = num2str(round([0,0.2,0.4,0.6,0.8,1]'*(max(mod_depth)-min(mod_depth)) + min(mod_depth)));
    b=colorbar('TickLabels',tick_labels);
    b.Label.String = 'Move - center hold (Hz)';
    colormap(inferno_list(1:end-19,:));
    ax = gca;
    plot([0,max([ax.XLim,ax.YLim])],[0,max([ax.XLim,ax.YLim])],'k--')
    formatForLee(gcf);
    xlabel('Center-hold evoked response (Hz)');
    ylabel('Move evoked response (Hz)');
    set(gca,'fontsize',14)
    
    figure(); hold on;
    plot(mod_depth,move_fr_center_fr_diff,'.');
    ax = gca;
    plot([max([ax.XLim,ax.YLim]),-max([ax.XLim,ax.YLim])],[-max([ax.XLim,ax.YLim]),max([ax.XLim,ax.YLim])],'k--')
    formatForLee(gcf);
    xlabel('Move - center hold no stim (Hz)');
    ylabel('Move - center hold stim (Hz)');
    set(gca,'fontsize',14)
%     
    figure(); hold on;
    plot(move_fr_baseline,move_fr,'.');
    ax = gca;
    plot([0,max([ax.XLim,ax.YLim])],[0,max([ax.XLim,ax.YLim])],'k--')
    formatForLee(gcf);
    xlabel('movement baseline response (Hz)');
    ylabel('stim during movement response (Hz)');
    set(gca,'fontsize',14)
    
%     errorbar(center_hold_fr(~stim_chan_mask),move_fr(~stim_chan_mask,1),...
%         move_std(~stim_chan_mask,1),move_std(~stim_chan_mask,1),...
%         center_hold_std(~stim_chan_mask),center_hold_std(~stim_chan_mask),...
%         '.','color',getColorFromList(1,0),'markersize',12);
    
%     errorbar(center_hold_fr(stim_chan_mask==0),move_fr(stim_chan_mask==0,2),...
%         move_std(stim_chan_mask==0,2),move_std(stim_chan_mask==0,2),...
%         center_hold_std(stim_chan_mask==0),center_hold_std(stim_chan_mask==0),...
%         '.','color',getColorFromList(1,1),'markersize',12);
    
%     errorbar(center_hold_fr(stim_chan_mask==1),move_fr(stim_chan_mask==1,1),...
%         move_std(stim_chan_mask==1,1),move_std(stim_chan_mask==1,1),...
%         center_hold_std(stim_chan_mask==1),center_hold_std(stim_chan_mask==1),...
%         '*','color',getColorFromList(1,0),'markersize',8);
    
%     errorbar(center_hold_fr(stim_chan_mask==1),move_fr(stim_chan_mask==1,2),...
%         move_std(stim_chan_mask==1,2),move_std(stim_chan_mask==1,2),...
%         center_hold_std(stim_chan_mask==1),center_hold_std(stim_chan_mask==1),...
%         '*','color',getColorFromList(1,1),'markersize',8);
%     


%% dPCA
    remove_stim_chan = 0; stim_chan = 25;
    bins_around_movement = [-5,12];
    
    normalize_mag = 0;
    plot_stim = 0;
    do_smooth = 0;
    
    color_list = getColorFromList(1,0:3);
    td_trim = td_all;
    if(remove_stim_chan==1)
        bad_units = find(td_trim(1).([array_name,'_unit_guide'])(:,1) == stim_chan);
        for trial = 1:length(td_trim)
            temp = td_trim(trial).([array_name '_spikes']);
            temp(:,bad_units) = [];
            td_trim(trial).([array_name '_spikes']) = temp;
            temp = td_trim(trial).([array_name '_unit_guide']);
            temp(bad_units,:) = [];
            td_trim(trial).([array_name '_unit_guide']) = temp;
        end
    end
    
    if(do_smooth)
        td_trim = smoothSignals(td_trim,struct('signals','LeftS1_spikes'));
    end
    
    td_trim = trimTD(td_trim,{'idx_goCueTime',bins_around_movement(1)},{'idx_goCueTime',bins_around_movement(2)});
    
    num_lines = size(td_trim(1).LeftS1_spikes,1);
    
    unique_tgt_dir = unique([td_trim.target_direction]);
   
    
    % split movement only and move-stim trials
    td_move = td_trim(isnan([td_trim.idx_stimTime]));
    td_stim = td_trim(~isnan([td_trim.idx_stimTime]));
    td_stim = td_stim([td_stim.idx_stimTime] > [td_stim.idx_goCueTime]);
    
    go_to_stim_idx = ceil(mean([td_stim.idx_stimTime] - [td_stim.idx_goCueTime]));
    stim_idx = abs(bins_around_movement(1)) + 1 + go_to_stim_idx+[0,1,2]; % bins with stims in them
    
     pca_params = []; pca_params.signals = {'LeftS1_spikes'}; pca_params.do_plot = 1;
    pca_params.marg_names = {'time','target'}; 
    [td_move,dpca_info] = runDPCA(td_move,'target_direction', pca_params);
    
    idx_plot = find(dpca_info.which_marg == 2,2,'first');
    
    td_move_avg = trialAverage(td_move,{'target_direction'});
    td_stim_avg = trialAverage(td_stim,{'target_direction'});

    if(~plot_stim) figure(); hold on; ax = gca; end
    % plot averages
    for i_trial = 1:numel(td_move_avg)
        if(plot_stim) figure('Position',[680 558 864 420]); subplot(1,2,1); hold on; ax = gca; end
        
        color_idx = find(td_move_avg(i_trial).target_direction == unique_tgt_dir);
        td_move_avg(i_trial).LeftS1_pca = td_move_avg(i_trial).LeftS1_spikes*dpca_info.W;
        td_stim_avg(i_trial).LeftS1_pca = td_stim_avg(i_trial).LeftS1_spikes*dpca_info.W;
        
        
        color_grad = [linspace(0,color_list(color_idx,1),num_lines)',linspace(0,color_list(color_idx,2),num_lines)',linspace(0,color_list(color_idx,3),num_lines)'];
        ax.ColorOrder = color_grad;
        ax.ColorOrderIndex = 1;
        x = [td_move_avg(i_trial).LeftS1_pca(1:end-1,idx_plot(1)), td_move_avg(i_trial).LeftS1_pca(2:end,idx_plot(1))];
        y = [td_move_avg(i_trial).LeftS1_pca(1:end-1,idx_plot(2)), td_move_avg(i_trial).LeftS1_pca(2:end,idx_plot(2))];
        plot(x',y','linewidth',2)
        plot(td_move_avg(i_trial).LeftS1_pca(1,idx_plot(1)),td_move_avg(i_trial).LeftS1_pca(1,idx_plot(2)),'.',...
            'color',color_grad(1,:),'markersize',35);
        
        formatForLee(gcf);
        xlabel(['PC',num2str(idx_plot(1))]);
        ylabel(['PC',num2str(idx_plot(2))]);
                  
       
        % plot stim average
        if(plot_stim)
            ax.ColorOrderIndex = 1;
            x = [td_stim_avg(i_trial).LeftS1_pca(1:end-1,idx_plot(1)), td_stim_avg(i_trial).LeftS1_pca(2:end,idx_plot(1))];
            y = [td_stim_avg(i_trial).LeftS1_pca(1:end-1,idx_plot(2)), td_stim_avg(i_trial).LeftS1_pca(2:end,idx_plot(2))];
            plot(x',y','linewidth',2,'linestyle','--');
            plot(td_stim_avg(i_trial).LeftS1_pca(1,idx_plot(1)),td_stim_avg(i_trial).LeftS1_pca(1,idx_plot(2)),...
                    '.','color',color_grad(1,:),'markersize',35); 
            % plot stim bins
            ax.ColorOrderIndex = stim_idx(1); 
            plot(td_stim_avg(i_trial).LeftS1_pca(stim_idx,idx_plot(1)),td_stim_avg(i_trial).LeftS1_pca(stim_idx,idx_plot(2)),...
                '*','markersize',12);          
            
            ax.ColorOrderIndex = stim_idx(1);
            plot(td_move_avg(i_trial).LeftS1_pca(stim_idx,idx_plot(1)),td_move_avg(i_trial).LeftS1_pca(stim_idx,idx_plot(2)),...
                '*','markersize',12); 
            % plot magnitude of projection onto each PC
            subplot(1,2,2); hold on
            if(normalize_mag)
                norm_factor = mean(abs(dpca_info.W));
            else
                norm_factor = ones(1,size(td_move_avg(i_trial).LeftS1_pca,2));
            end
            data_1 = td_move_avg(i_trial).LeftS1_pca(stim_idx(1)-1,:);
            data_2 = td_move_avg(i_trial).LeftS1_pca(stim_idx(end),:);
            move_magnitude = abs(data_1-data_2);
            move_magnitude = move_magnitude./norm_factor; 

            data_1 = td_stim_avg(i_trial).LeftS1_pca(stim_idx(1)-1,:);
            data_2 = td_stim_avg(i_trial).LeftS1_pca(stim_idx(end),:);
            stim_magnitude = abs(data_1-data_2);
            stim_magnitude = stim_magnitude./norm_factor; 

            plot(move_magnitude,'color',color_list(i_trial,:),'linewidth',2)
            plot(stim_magnitude,'--','color',color_list(i_trial,:),'linewidth',2)
            plot(idx_plot,stim_magnitude(idx_plot),'.','color',color_list(i_trial,:),'markersize',20)
            ylabel('Magnitude'); xlabel('PC');
            ylim([0,3.5])
        end
    end
    
    % get difference between stim trials and trial averaged path for each
    % bin
    %%
    subplot_idx = [6,2,4,8];
    offset = [-1,1]/500;
    color_list = {'k','r'};
    plot_magnitude = 1; % else plot non-absolute value distance
    
    pos_plot_idx = stim_idx(3);
    
%     f_dist = figure();
    f_pos = figure();
%     stim_pos_data = nan(numel(td_stim),numel(stim_idx)+1,size(td_move_avg(1).LeftS1_pca,2));
%     move_pos_data = nan(numel(td_move),numel(stim_idx)+1,size(td_move_avg(1).LeftS1_pca,2));
    for i = 1:2 % td_move then td_stim
        if(i==1) td_test = td_move; end
        if(i==2) td_test = td_stim; end
        
        tgt_idx = nan(numel(td_test),1);
        distance_data = nan(numel(td_test),size(td_move_avg(1).LeftS1_pca,1),2); % in tgt_dir space, and in all spaces
        dimension_distance_data = nan(numel(td_test),size(td_move_avg(1).LeftS1_pca,1),size(td_move_avg(1).LeftS1_pca,2));
        pos_data = nan(numel(td_test),numel(pos_plot_idx),size(td_move_avg(1).LeftS1_pca,2));
        
        for i_trial = 1:numel(td_test)
            % get tgt_idx
            tgt_idx(i_trial) = find(td_test(i_trial).target_direction == unique_tgt_dir);

            % project trial into dpca space
            td_test(i_trial).LeftS1_pca = td_test(i_trial).LeftS1_spikes*dpca_info.W;

            % compute distance in target direction space (idx_plot) and in all
            % spaces
            dimension_distance_data(i_trial,:,:) = (td_test(i_trial).LeftS1_pca - td_move_avg(tgt_idx(i_trial)).LeftS1_pca);
            distance_data(i_trial,:,1) = sqrt(sum((td_test(i_trial).LeftS1_pca(:,idx_plot) - td_move_avg(tgt_idx(i_trial)).LeftS1_pca(:,idx_plot)).^2,2));
            distance_data(i_trial,:,2) = sqrt(sum((td_test(i_trial).LeftS1_pca(:,:) - td_move_avg(tgt_idx(i_trial)).LeftS1_pca(:,:)).^2,2));
            pos_data(i_trial,:,:) = td_test(i_trial).LeftS1_pca(pos_plot_idx,:);
        end

        for i_tgt = 1:numel(unique_tgt_dir)
%             figure(f_dist);
%             subplot(3,3,subplot_idx(i_tgt)); hold on
%             x = ((1:1:size(distance_data,2))-1+bins_around_movement(1))*td_test(1).bin_size;
%             y = mean(distance_data(tgt_idx==i_tgt,:,1),1);
%             err = std(distance_data(tgt_idx==i_tgt,:,1),1);
%             errorbar(x+offset(i),y,err,'color',color_list{i});
%             if(i==2) % stim trials)
%                 plot(x(stim_idx),y(stim_idx),'*','color',color_list{i},'markersize',8)
%             end
%             
%             xlabel('Time relative to movement onset (s)');
%             ylabel('Distance from move average');
%             xlim([x(1),x(end)]);
            
            figure(f_pos);
            subplot(3,3,subplot_idx(i_tgt)); hold on
            plot(pos_data(tgt_idx==i_tgt,:,idx_plot(1)),pos_data(tgt_idx==i_tgt,:,idx_plot(2)),'.','color',color_list{i})
            plot(mean(pos_data(tgt_idx==i_tgt,:,idx_plot(1))),mean(pos_data(tgt_idx==i_tgt,:,idx_plot(2))),'.','color',color_list{i},'markersize',30)
        end
    end