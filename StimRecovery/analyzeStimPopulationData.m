%% set initial parameters

    input_data.folderpath = 'D:\Lab\Data\StimDuringMove\Han_20191004_CObump_stimDuringTask\';

%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\6250-002085\SN 6250-002085.cmp';
    
    input_data.date = '20191004';

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
    
    
    array_name = input_data.array(6:end);
%% make cds and trial data
    cd(input_data.folderpath)

    file_name = dir('*Extracted.nev*');
    
    td_all = [];
    bin_size = 0.05; % s
    
    for i_file = 1:numel(file_name)
        cds = commonDataStructure();
        cds.file2cds(strcat(input_data.folderpath,file_name(i_file).name),input_data.array,input_data.monkey,input_data.ranBy,...
            input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
        
        
        params.event_list = {'goCueTime';'tgtDir'};
        params.trial_results = {'R','F'};
        params.extra_time = [1,4];
        params.include_ts = 0;
        params.exclude_units = [0,255];
        
        td_temp = parseFileByTrial(cds,params);
        td_temp = getNorm(td_temp,'vel');
        td_temp = getMoveOnsetAndPeak(td_temp);
        td_temp = removeBadTrials(td_temp);

        % combine trial datas -- deal with spike data that is different
        % sizes (neuron dropped out or something weird).
        td_all = combineTDs(td_all,td_temp);
    end

    td_all = binTD(td_all,ceil(bin_size/td_all(1).bin_size));
        
    
%% plot neural activity during movements to different targets (with and without stim)
	do_plot = 1;
    
%     subplot_idx = [7,8,9,6,3,2,1,4];
    subplot_idx = [6,2,4,8];
    td_move = trimTD(td_all,{'idx_tgtOnTime',ceil(-0.3/bin_size)},{'idx_tgtOnTime',ceil(0.6/bin_size)});
    tgt_dirs = unique([td_move.target_direction]);
    if(do_plot)
        for i_unit = 1:size(td_all(1).([array_name,'_spikes']),2)
            figure();
            for i_tgt = 1:numel(tgt_dirs)
                % get average PSTH then plot
                trial_mask = [td_move.target_direction] == tgt_dirs(i_tgt);
                td_tgt = td_move(trial_mask == 1);
                % get average (trialAverage fails for older versions of matlab)
                spike_count = zeros(size(td_tgt(1).([array_name,'_spikes']),1),1);
                for i_trial = 1:numel(td_tgt)
                    spike_count = spike_count + td_tgt(1).([array_name,'_spikes'])(:,i_unit);
                end
                spike_count = spike_count./numel(td_tgt);

                subplot(3,3,subplot_idx(i_tgt))
                plot(spike_count/td_move(1).bin_size);
            end
        end
    end
    
td_move_all = td_move;

   
%% plot emg activity during movements to a given target for example trials
    td_move = trimTD(td_all,{'idx_movement_on',ceil(-1/bin_size)},{'idx_movement_on',ceil(3.0/bin_size)});

    tgt_dirs = unique([td_move.target_direction]);
    
    for i_tgt = 1:numel(tgt_dirs)
        f=figure('Position',[2157,108,560,861]);
        td_move_tgt = td_move([td_move.target_direction] == tgt_dirs(i_tgt));
        td_move_tgt = trialAverage(td_move_tgt);
        trial_idx = datasample(1:1:numel(td_move_tgt),1);
        x_data = ((1:1:size(td_move_tgt(trial_idx).emg,1))-td_move_tgt(trial_idx).idx_movement_on)*td_move_tgt(trial_idx).bin_size;
        for i_emg = 1:numel(td_move(i_tgt).emg_names)
            
            subplot(numel(td_move(i_tgt).emg_names)/2,2,i_emg);
            
            plot(x_data,td_move_tgt(trial_idx).emg(:,i_emg),'k','linewidth',2)
            ylabel(td_move(i_tgt).emg_names{i_emg}(5:end));
            set(gca,'fontsize',12)
            formatForLee(gcf)
            if(i_emg == numel(td_move(i_tgt).emg_names))
                xlabel('Time after movement on (s)');
            end
            xlim([x_data(1),x_data(end)])
        end
        
        f.Name = ['Pop_20201023_EMG_tgtDir',num2str(tgt_dirs(i_tgt))];
        saveFiguresLIB(f,input_data.folderpath,f.Name);
    end
%% define recovery metric and run on a bunch of array data's
    model_name = 'dpca';
    sqrt_transform = 0;
    
    % arrayData files need to be in same folderpath as task related data
    cd(input_data.folderpath)
    array_data_files = dir('*arrayData*');
    
    % for each file, load array data, trim TD appropriately, compute dpc's,
    % then get metric and store
    td_move_all = trimTD(td_all,{'idx_movement_on',ceil(-0.5/bin_size)},{'idx_movement_on',ceil(0.7/bin_size)});
    dist_to_baseline = [];
    recov_ts = [];
    dir_proj = [];
    
    f_stim_mag = figure();hold on;
    for i_file = 1:numel(array_data_files)
        load([input_data.folderpath,array_data_files(i_file).name]);
        td_move = td_move_all;
        if(sqrt_transform)
            td_move = sqrtTransform(td_move);
        end
        
        [td_move,arrayData] = matchChannelListsTDArrayData(td_move,arrayData,array_name);
        
        if(numel(arrayData) > 20)
            % dpca or pca based on model_name
            % function handles doing pca or dpca with corresponding parameters
            [td_move_pca, model_info,idx_proj] = dimReduceTDWrapper(td_move, model_name, array_name);

            % get neural state for each dimension for each stim
            [neural_proj, arrayData_rebin] = getArrayDataProjection(arrayData, model_info, bin_size*1000, sqrt_transform,0);
            mean_neural_proj = squeeze(mean(neural_proj,2));
            
            % get recovery metric and store
            recovery_data = getPopulationRecovery(neural_proj, arrayData_rebin{1}.binEdges{1});

            recov_ts(i_file) = recovery_data.all_dim - 100; % 100 ms stim length
            recov_tau(i_file) = 1/recovery_data.all_dim_fit.b;
            dist_to_baseline(i_file,:) = recovery_data.dist_to_baseline;

            % get projection length for each dimension
            baseline_idx = find(arrayData_rebin{1}.binEdges{1} < 0);
            baseline_idx(end-1:end) = [];
            stim_idx = find(arrayData_rebin{1}.binEdges{1} >= 0);
            stim_idx = stim_idx(1:3);
            mag_data = (mean_neural_proj-mean(mean_neural_proj(:,baseline_idx),2))./(sqrt(sum(model_info.W.^2))');

            baseline_mag = mean(mag_data(:,baseline_idx),2);

            stim_mag = max(mag_data-baseline_mag,[],2)./(sqrt(sum(model_info.W.^2))');

            % get projection direction            
            dir_proj(i_file,:) = atan2(mean_neural_proj(idx_proj(2),stim_idx)-mean(mean_neural_proj(idx_proj(2),baseline_idx)),...
                mean_neural_proj(idx_proj(1),stim_idx)-mean(mean_neural_proj(idx_proj(1),baseline_idx)));
            mag_proj(i_file,:) = sqrt(sum((mean_neural_proj(idx_proj,stim_idx)-mean(mean_neural_proj(idx_proj,baseline_idx),2)).^2));
            
            % plot stim_mag across dims for this stimulation
            figure(f_stim_mag);
            plot(stim_mag);

            % plot distance from baseline
            f=figure(); hold on;
            f.Name = [array_data_files(i_file).name(1:end-4),'_distBaseline'];
            bin_edges = arrayData_rebin{1}.binEdges{1};

            plot(bin_edges(1:end-1) + mode(diff(bin_edges))/2,dist_to_baseline(i_file,:));
            plot([0,0],[0,5],'k--','linewidth',2)
            plot([100,100],[0,5],'k--','linewidth',2)
            plot([recov_ts(i_file)+100,recov_ts(i_file)+100],[0,0],'r.','markersize',20)
            
            
            saveFiguresLIB(f,input_data.folderpath, f.Name);
        end
            
    end
    
    
%% do dPCA and visualize
    model_name = 'dpca';
    sqrt_transform = 0;
    
    plot_stim = 0;
    wave_idx = 1;
   
    td_move_all = trimTD(td_all,{'idx_movement_on',ceil(-0.3/bin_size)},{'idx_movement_on',ceil(0.5/bin_size)});
    td_move = td_move_all;
    
    if(sqrt_transform == 1)
       td_move = sqrtTransform(td_move); 
    end
    
    [td_move,arrayData] = matchChannelListsTDArrayData(td_move,arrayData,array_name);
    
    color_list = inferno(9);
    
    unique_tgt_dir = unique([td_move.target_direction]);
    % function handles doing pca or dpca with corresponding parameters
    [td_move_pca,model_info,idx_plot] = dimReduceTDWrapper(td_move, model_name, array_name);
        
    td_move_pca_avg = trialAverage(td_move_pca,{'target_direction'});

    move_magnitude = [];
    f_proj = figure(); hold on;
    if(plot_stim) subplot(1,2,1); hold on; f_proj.Position = [694 543 1069 430]; end 
    for i_trial = 1:numel(td_move_pca_avg)
        color_idx = find(td_move_pca_avg(i_trial).target_direction == unique_tgt_dir);
        td_move_pca_avg(i_trial).LeftM1_pca = (td_move_pca_avg(i_trial).([array_name,'_spikes'])-model_info.mu)*model_info.W;
        
        plot(td_move_pca_avg(i_trial).LeftM1_pca(:,idx_plot(1)),td_move_pca_avg(i_trial).LeftM1_pca(:,idx_plot(2)),'color',color_list(color_idx,:),'linewidth',2);
        plot(td_move_pca_avg(i_trial).LeftM1_pca(1,idx_plot(1)),td_move_pca_avg(i_trial).LeftM1_pca(1,idx_plot(2)),'.','color',color_list(color_idx,:),'markersize',20);

        td_move_dir = td_move([td_move.target_direction] == unique_tgt_dir(i_trial));
        for i_single = 1:numel(td_move_dir)
            pca_temp = (td_move_dir(i_single).([array_name,'_spikes'])-model_info.mu)*model_info.W;
            plot(pca_temp(end,idx_plot(1)),pca_temp(end,idx_plot(2)),'.','color',color_list(color_idx,:))
        end
        
        formatForLee(gcf);
        xlabel(['dPC',num2str(idx_plot(1))]);
        ylabel(['dPC',num2str(idx_plot(2))]);
        
        data_1 = td_move_pca_avg(i_trial).LeftM1_pca(7,:);
        data_2 = td_move_pca_avg(i_trial).LeftM1_pca(14,:);
        move_magnitude(i_trial,:) = abs(data_1-data_2);
    end
    
    max_stim_time = 0.8;
    
    if(plot_stim && exist('arrayData') > 0) 
        bin_edges = -bin_size:bin_size:max_stim_time;
        num_bins = numel(bin_edges) - 1; % stim time of stim data to project
        spike_matrix = zeros(num_bins,numel(arrayData));
        mag_data = [];
        for arr_idx = 1:numel(arrayData)
            chan(arr_idx) = arrayData{arr_idx}.CHAN_REC;
            spike_matrix(:,arr_idx) = histcounts(arrayData{arr_idx}.spikeTrialTimes{wave_idx},bin_edges)/arrayData{arr_idx}.numStims(wave_idx);

        end
        % sqrt transform if requested
        if(sqrt_transform)
            spike_matrix = sqrt(spike_matrix);
        end
        
        % project into pca space found above
        pca_spike_matrix = (spike_matrix-model_info.mu)*model_info.W;
        
        plot(pca_spike_matrix(:,idx_plot(1)), pca_spike_matrix(:,idx_plot(2)),'color',getColorFromList(1,1),'linewidth',3,'linestyle','-','marker','.','markersize',20)
        plot(pca_spike_matrix(1,idx_plot(1)), pca_spike_matrix(1,idx_plot(2)),'s','markersize',12,'linewidth',2,'color',getColorFromList(1,1))

        data_1 = pca_spike_matrix(1,:);
        data_2 = pca_spike_matrix(find(bin_edges >= 0.1,1,'first')-1,:);
        mag_data = abs(data_1-data_2)./(sqrt(sum(model_info.W.^2)));
        
        subplot(1,2,2); hold on;
        plot(mean(move_magnitude)','k','linewidth',1.5);
        hold on
        plot(mag_data,'--','color',getColorFromList(1,1),'linewidth',2);
        plot(idx_plot,mag_data(idx_plot),'.','color',getColorFromList(1,1),'markersize',20)
        l=legend('movement','stim'); set(l,'box','off');
        formatForLee(gcf)
        xlabel('dPC'); ylabel('Magnitude');
    end


%% null vs potent space analysis -- where does stim end up?
    
    model_name = 'pca';
    sqrt_transform = 0;
    
    td_move_all = trimTD(td_all,{'idx_movement_on',ceil(-1/bin_size)},{'idx_movement_on',ceil(1.5/bin_size)});
    td_move = td_move_all;
    
    if(sqrt_transform == 1)
       td_move = sqrtTransform(td_move); 
    end
    
    [td_move,arrayData] = matchChannelListsTDArrayData(td_move,arrayData,array_name);

    pot_params = []; pot_params.in_signals = [array_name,'_spikes'];
    pot_params.out_signals = ['emg'];
    pot_params.in_dims = 16; pot_params.out_dims = 6;
    
    
    [td_move_pot,pca_info] = getPotentSpace(td_move,pot_params);
    
    
    % ok so plot mean null and potent activity against time for each tgt
    % null activity should increase prior to potent activity? it isn't a
    % delay period task though for the record
    td_move_pot_avg = trialAverage(td_move_pot,{'target_direction'});
    unique_tgt_dir = unique([td_move_pot_avg.target_direction]);
    
    figure();
    subplot_idx = [7,8,9,6,3,2,1,4];
    for i_tgt = 1:numel(unique_tgt_dir)
        subplot(3,3,subplot_idx(i_tgt))
        % get neural distance from baseline for potent and null space
        baseline_idx = 1:15;
        baseline_pot = mean(td_move_pot_avg(i_tgt).LeftM1emg_potent(baseline_idx,:));
        baseline_null = mean(td_move_pot_avg(i_tgt).LeftM1emg_null(baseline_idx,:));
        
        dist_pot = td_move_pot_avg(i_tgt).LeftM1emg_potent - baseline_pot;
        dist_null = td_move_pot_avg(i_tgt).LeftM1emg_null - baseline_null;
        
        x_data = ((1:1:size(td_move_pot_avg(i_tgt).LeftM1emg_null,1))-td_move_pot_avg(i_tgt).idx_movement_on)*td_move_pot_avg(i_tgt).bin_size;
        
        plot(x_data,sqrt(sum(dist_pot.^2,2)),'color',getColorFromList(1,0));
        hold on;
        plot(x_data,sqrt(sum(dist_null.^2,2)),'color',getColorFromList(1,1))
        
        formatForLee(gcf);
        ylabel('Euclidean dist from baseline');
        xlabel('Time from movement on (s)');
        
        if(i_tgt == 7)
           l=legend('Potent','Null');
        end
    end
    
    % ok so let's project stim data into the null and potent space, then
    % quantify how far the neural state moves in each dimension
    [neural_proj, arrayData_rebin] = getArrayDataProjection(arrayData, pca_info, bin_size*1000, sqrt_transform,1);
    
    neural_proj.mean_pot = squeeze(mean(neural_proj.pot_proj,2));
    neural_proj.mean_null = squeeze(mean(neural_proj.null_proj,2));
    
    baseline_pot = mean(neural_proj.mean_pot(:,1:15),2);
    baseline_null = mean(neural_proj.mean_null(:,1:15),2);
    
    pot_dist = sqrt(sum((neural_proj.mean_pot - baseline_pot).^2));
    null_dist = sqrt(sum((neural_proj.mean_null -baseline_null).^2));
    
    figure();
    x_data = arrayData_rebin{1}.binEdges{1}(1:end-1) + mode(diff(arrayData_rebin{1}.binEdges{1}));
    plot(x_data,pot_dist,'color',getColorFromList(1,0)); hold on;
    plot(x_data,null_dist,'color',getColorFromList(1,1));
    
    xlabel('Time from stim onset (s)');
    ylabel('Euclidean distance from baseline');
    formatForLee(gcf);
    set(gca,'fontsize',14)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    