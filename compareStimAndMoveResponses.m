%% set initial parameters

    input_data.folderpath = 'D:\Lab\Data\StimPDs\Han_20190927\';

%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\6250-002085\SN 6250-002085.cmp';
    
    input_data.date = '20190923';

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
    
    cds = commonDataStructure();
    cds.file2cds(strcat(input_data.folderpath,file_name(1).name),input_data.array,input_data.monkey,input_data.ranBy,...
        input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
    
    
%% make trial data
    BIN_SIZE = 0.02; NUM_BINS = 5;
    
    params.event_list = {'goCueTime'};
    params.trial_results = {'R'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [255];

    move_onset_params.pre_move_thresh = 1000;
    move_onset_params.min_s = 3;
    move_onset_params.start_idx_offset = -10;
    move_onset_params.max_rt_offset = 400;
    
    td_all = parseFileByTrial(cds,params);
    td_all = stripSpikeSorting(td_all);
    td_all = getSpeed(td_all);
%     td_all = getMoveOnset(td_all,move_onset_params);
    td_all = removeBadTrials(td_all);

%     if(td_all(1).bin_size < 1)
        % set it to 50ms
        td_all = binTD(td_all,ceil(BIN_SIZE/td_all(1).bin_size));
%     end
    
%% get correlation during movement epoch
    corr_params = [];
    corr_params.signals = {'LeftS1_spikes'};
    td_move = trimTD(td_all,{'idx_goCueTime',-5},{'idx_goCueTime',15});
    td_move = softNormalize(td_move);

%     [td_move_avg,cond_idx] = trialAverage(td_move,{'target_direction'});
%     td_move_mean_sub = subtractConditionMean(td_move);
    
    [rho_move,sort_idx] = pairwiseCorr(td_move,corr_params);
    
    
%% load in array data
    

%% get correlation matrix for stim data
    
    stim_bin_size = BIN_SIZE*1000; % ms
    wave_idx = 4;
    % for each stim trial, bin data according to bin size above and store
    % making a # neurons x # bins x # stim trials matrix
    num_neurons = numel(arrayData);
    num_stims = arrayData{1}.numStims(wave_idx);
    
    num_bins_before_0 = floor(abs(min(arrayData{1}.binEdges{1}))/(stim_bin_size));
    num_bins_after_0 = floor(abs(max(arrayData{1}.binEdges{1}))/(stim_bin_size));
    num_bins = num_bins_before_0 + num_bins_after_0;
    bin_edges = (-num_bins_before_0:1:num_bins_after_0)*stim_bin_size;
    
    spike_matrix = zeros(num_neurons,num_bins,num_stims);
    chan = zeros(num_neurons,1);
    wave_sent_list = arrayData{1}.WAVEFORM_SENT(1:100:end);
    for arr_idx = 1:numel(arrayData)
        chan(arr_idx) = arrayData{arr_idx}.CHAN_REC;
        for stim_idx = 1:num_stims
            spike_mask = arrayData{arr_idx}.stimData{wave_idx} == stim_idx;
            spike_matrix(arr_idx,:,stim_idx) = histcounts(arrayData{arr_idx}.spikeTrialTimes{wave_idx}(spike_mask),bin_edges/1000);
        end
    end
    
    % compute pairwise correlation and store as rho_stim
    rho_stim = zeros(num_neurons,num_neurons,3); % pre stim, during stim, post stim
    spikes_pre_stim = spike_matrix(:,num_bins_before_0-NUM_BINS:num_bins_before_0-1,:);
    rho_stim(:,:,1) = corr(reshape(spikes_pre_stim,num_neurons,NUM_BINS*num_stims)') + -1*eye(num_neurons,num_neurons);
    spikes_dur_stim = spike_matrix(:,num_bins_before_0+1:num_bins_before_0+NUM_BINS,:);
    rho_stim(:,:,2) = corr(reshape(spikes_dur_stim,num_neurons,NUM_BINS*num_stims)') + -1*eye(num_neurons,num_neurons);
    stim_dur = ceil(0.1/BIN_SIZE);
    spikes_post_stim = spike_matrix(:,num_bins_before_0+2+stim_dur:num_bins_before_0+2+stim_dur+NUM_BINS-1,:);
    rho_stim(:,:,3) = corr(reshape(spikes_post_stim,num_neurons,NUM_BINS*num_stims)') + -1*eye(num_neurons,num_neurons);
    
    % remove neurons from rho_move that are on different electrodes
    
    rho_move_trim = rho_move(chan,:);
    rho_move_trim = rho_move_trim(:,chan);
%     rho_move_trim = rho_move;
    
    figure('Position',[2037 248 1748 262])
    subplot(1,4,1)
    
    imagesc(rho_move_trim); colorbar
    title('move')
    subplot(1,4,2)
    
    imagesc(rho_stim(:,:,1)); colorbar
    title('prestim')
    subplot(1,4,3)
    
    imagesc(rho_stim(:,:,2)); colorbar
    title('stim')
    subplot(1,4,4)
    
    imagesc((rho_stim(:,:,2)-rho_stim(:,:,1))); colorbar
    title('stim - prestim')
    
    
%% trim td_move channels based on those in arrayData
    chan = [];
    for arr_idx = 1:numel(arrayData)
        chan(arr_idx) = arrayData{arr_idx}.CHAN_REC;
    end

    for i_trial = 1:numel(td_move)
        chan_mask = any(td_move(i_trial).LeftS1_unit_guide(:,1) == chan,2);
        td_move(i_trial).LeftS1_spikes = td_move(i_trial).LeftS1_spikes(:,chan_mask==1);
        td_move(i_trial).LeftS1_unit_guide = td_move(i_trial).LeftS1_unit_guide(chan_mask==1,:);
    end


%% can we look at the effect of stim in a movement predictive space
% build pca space, project stim data into it
    plot_stim = 1;
    wave_idx = 3;
    idx_plot = [2,3];
    unique_tgt_dir = unique([td_move.target_direction]);
    color_list = inferno(9);
    pca_params = []; pca_params.signals = {'LeftS1_spikes'};
    
    [td_move_pca,pca_info] = getPCA(td_move, pca_params);
    td_move_pca = getPCA(td_move,pca_info);
    td_move_pca_avg =  trialAverage(td_move_pca,{'target_direction'});
    
    f=figure();  
    if(plot_stim) subplot(1,2,1); f.Position = [680 558 1018 420]; end 
    hold on;
    move_magnitude = zeros(numel(td_move_pca_avg),size(td_move_pca_avg(1).LeftS1_pca,2));
    for i_trial = 1:numel(td_move_pca_avg)
        color_idx = find(td_move_pca_avg(i_trial).target_direction == unique_tgt_dir);
        
        plot(td_move_pca_avg(i_trial).LeftS1_pca(:,idx_plot(1)),td_move_pca_avg(i_trial).LeftS1_pca(:,idx_plot(2)),'color',color_list(color_idx,:),'linewidth',2);
        plot(td_move_pca_avg(i_trial).LeftS1_pca(1,idx_plot(1)),td_move_pca_avg(i_trial).LeftS1_pca(1,idx_plot(2)),'.','color',color_list(color_idx,:),'markersize',20);
       
        formatForLee(gcf);
        xlabel(['PC',num2str(idx_plot(1))]);
        ylabel(['PC',num2str(idx_plot(2))]);
        
        data_1 = td_move_pca_avg(i_trial).LeftS1_pca(1,:);
        data_2 = td_move_pca_avg(i_trial).LeftS1_pca(end,:);
        move_magnitude(i_trial,:) = sqrt((data_1-data_2).^2 + (data_1-data_2).^2);
    end

%     for i_trial = 1:numel(td_move_pca)
%         color_idx = find(td_move_pca(i_trial).target_direction == unique_tgt_dir);
%         
%         plot(td_move_pca(i_trial).LeftS1_pca(end,2),td_move_pca(i_trial).LeftS1_pca(end,3),'.','color',color_list(color_idx,:),'markersize',10);
%         
%     end
    % project stim data into pca space
    max_stim_time = 0.10;
    
    num_bins = floor(max_stim_time/BIN_SIZE) + 1; % 200 ms of stim data to project
    bin_edges = -BIN_SIZE:BIN_SIZE:max_stim_time;
    spike_matrix = zeros(num_bins,numel(arrayData));
    stim_magnitude = [];
    if(plot_stim && exist('arrayData') > 0) 
        for arr_idx = 1:numel(arrayData)
            chan(arr_idx) = arrayData{arr_idx}.CHAN_REC;
            spike_matrix(:,arr_idx) = histcounts(arrayData{arr_idx}.spikeTrialTimes{wave_idx},bin_edges)/arrayData{arr_idx}.numStims(wave_idx);
            % soft normalize spikes
            spike_matrix(:,arr_idx) = spike_matrix(:,arr_idx)./(range(spike_matrix(:,arr_idx))+5);
        end
        
        % project into pca space found above
        pca_spike_matrix = (spike_matrix - repmat(pca_info.mu,size(spike_matrix,1),1))*pca_info.w;
        
        plot(pca_spike_matrix(:,idx_plot(1)), pca_spike_matrix(:,idx_plot(2)),'color',getColorFromList(1,1),'linewidth',3,'linestyle','--')
        plot(pca_spike_matrix(1,idx_plot(1)), pca_spike_matrix(1,idx_plot(2)),'.','markersize',20,'color',getColorFromList(1,1))

        data_1 = pca_spike_matrix(1,:);
        data_2 = pca_spike_matrix(end,:);
        stim_magnitude = sqrt((data_1-data_2).^2 + (data_1-data_2).^2);
        
        subplot(1,2,2)
        plot(mean(move_magnitude)','k','linewidth',1.5);
        hold on
        plot(stim_magnitude,'--','color',getColorFromList(1,1),'linewidth',2);
        plot(idx_plot,stim_magnitude(idx_plot),'.','color',getColorFromList(1,1),'markersize',20)
        l=legend('movement','stim'); set(l,'box','off');
        formatForLee(gcf)
        xlabel('PC'); ylabel('Magnitude');
    end


%% investigate effect of stimulating a single neuron (stimulated channel esque)
    % plot weights for a single neuron across dimensions
    unit = 1;
    figure()
    plot(pca_info.w(unit,:))









