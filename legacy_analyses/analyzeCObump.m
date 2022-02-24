%% set initial parameters

    input_data.folderpath = 'D:\Lab\Data\Kramer_RW\current\';

%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\6250-002085\SN 6250-002085.cmp';
    
    %Kramer Right S1
    mapFileName = 'R:\limblab\lab_folder\Lab-Wide Animal Info\Implants\Blackrock Array Info\Array Map Files\6251-0922\6251-0922.cmp';

    %Kramer Left S1
%     mapFileName = 'R:\limblab\lab_folder\Lab-Wide Animal Info\Implants\Blackrock Array Info\Array Map Files\1024-0589\1024-0589.cmp';
    
    input_data.date = '20130404';

    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyKramer';
    input_data.ranBy = 'ranByTucker';
    input_data.lab = 6;
    input_data.mapFile = strcat('mapFile',mapFileName);
    input_data.task = 'taskRW';

    pwd=cd;
    input_data.fileNum = 1;
  
    input_data.center_y = -33;
    input_data.center_x = 3;
    input_data.num_bins = 8;
    
%% make cds
    cd(input_data.folderpath)

    file_name = dir('*1.nev*');
    
    params.event_list = {'goCueTime';'stimTime'};
    params.trial_results = {'R'};
    params.extra_time = [1,2];
    params.include_ts = 1;
    params.exclude_units = [0,255];

    move_onset_params.pre_move_thresh = 1000;
    move_onset_params.min_s = 3;
    move_onset_params.start_idx_offset = -10;
    move_onset_params.max_rt_offset = 400;
    
    td_all = [];
    for i = 1%:numel(file_name)
        cds = commonDataStructure();
        cds.file2cds(strcat(input_data.folderpath,file_name(i).name),input_data.array,input_data.monkey,input_data.ranBy,...
            input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
       
%         td_all = parseFileByTrial(cds, params);
        
    end
    
%     td_all = removeBadTrials(td_all);
%     if(td_all(1).bin_size < 0.05)
%         % set it to 50ms
%         td_all = binTD(td_all,ceil(0.01/td_all(1).bin_size));
%     end
    

%% get correlation during movement epoch
    corr_params = [];
    corr_params.signals = {'LeftS1_spikes'};
    td_move = trimTD(td_all,{'idx_goCueTime',0},{'idx_goCueTime',10});
    td_move = softNormalize(td_move);

%     avg_params = []; avg_params.conditions = 'all';
%     [td_move_avg,cond_idx] = trialAverage(td_move,{'target_direction'});
    td_move_mean_sub = subtractConditionMean(td_move);
    
    [rho,sort_idx] = pairwiseCorr(td_move_mean_sub,corr_params);


    
%% get PDs
    
    pd_params = [];
    pd_params.out_signals = 'LeftS1_spikes';
    pd_params.bootForTuning = 1;
    pd_params.num_boots = 100;
    pd_params.move_corr = 'vel';
    
    pd_all = getTDPDs(td_all,pd_params);
    
    
    
%% get PDs for bump

    td_bump = trimTD(td_all,{'idx_bumpTime',0},{'idx_bumpTime',300});

    pd_params = [];
    pd_params.out_signals = 'LeftS1_spikes';
    pd_params.bootForTuning = 0;
    pd_params.num_boots = 1;
    pd_params.move_corr = 'vel';
    
    pd_bump = getTDPDs(td_bump,pd_params);
    
%% get PDs for move
    td_move = trimTD(td_all,{'idx_goCueTime',0},{'idx_endTime',0});

    pd_params = [];
    pd_params.out_signals = 'LeftS1_spikes';
    pd_params.bootForTuning = 0;
    pd_params.num_boots = 1;
    pd_params.move_corr = 'vel';
    
    pd_move = getTDPDs(td_move,pd_params);

      
%% visualize
%% heatmap for preferred directions
disp('start')
    mapData = loadMapFile(mapFileName);
    
    optsPD.MAKE_BAR_PLOT = 1;

    optsPD.PLOT_CHANNELS = [1:96];
    optsPD.STIM_CHANNEL = 41;

    optsPD.MAX_RATIO = 1;
    optsPD.MIN_RATIO = -1;
    optsPD.LOG_SCALE = 0;
    optsPD.LOG_PARAM = 9;

    optsPD.FIGURE_SAVE = 0;
    optsPD.FIGURE_DIR = input_data.folderpath;
    optsPD.FIGURE_PREFIX = 'Han_20190924';
    
    [heatmapPD] = plotHeatmapsPD(td_all,pd_all,mapData,optsPD);

    
%% raster plots for each direction for a neuron
    td = td_all(isnan([td_all.idx_bumpTime]) & [td_all.result] == 'R');
    window_plot = [-0.5,0.5]; % s 
    window_pd = [0.,0.125]; % s
    speed_data = []; mean_speed_data = [];
    speed_x_data = window_plot(1):td(1).bin_size:window_plot(2);
    ax_list = [];
% for unit_idx = 1:72
    unit_idx = 40; % 1,5,32,33,39,40,58,62 looks good for Duncan's file
    try
        optsPlot = []; optsSave = [];
        optsPlot.MAKE_FIGURE = 0;
        optsPlot.MARKER_STYLE = 'line';
        optsPlot.LINE_WIDTH = 0.8;
        
        unique_tgt_dirs = unique([td.target_direction]);
        tgt_dir_mean_fr = [];
        tgt_dir_std_fr = [];
        trial_fr_all = [];
        trial_tgt_dir_all = [];
        subplot_map = [6,3,2,1,4,7,8,9];
        f_raster = figure('Position',[680 237 1081 741]);
        f_speed = figure('Position',[680 237 1081 741]);
        for i_tgt = 1:numel(unique_tgt_dirs)
            
            x_data = [];
            y_data = [];
            trial_fr = [];
            td_dir = td([td.target_direction] == unique_tgt_dirs(i_tgt));

            for i_trial = 1:numel(td_dir)
                % td_dir spike times are relative to trial onset
                spike_times_movement_on = td_dir(i_trial).LeftS1_ts{unit_idx} - ...
                    (td_dir(i_trial).idx_movement_on-td_dir(i_trial).idx_startTime)*td_dir(i_trial).bin_size;

                spike_mask = spike_times_movement_on >= window_plot(1) & spike_times_movement_on <= window_plot(2);

                x_data(end+1:end+sum(spike_mask)) = spike_times_movement_on(spike_mask);
                y_data(end+1:end+sum(spike_mask)) = i_trial;
                
                fr_mask = spike_times_movement_on >= window_pd(1) & spike_times_movement_on <= window_pd(2);
                trial_fr(i_trial) = sum(fr_mask)/diff(window_pd);
                trial_fr_all(end+1) = trial_fr(i_trial);
                trial_tgt_dir_all(end+1) = unique_tgt_dirs(i_tgt);
                
                % also want to collect speed data
                idx_speed = td_dir(i_trial).idx_movement_on + ceil(window_plot/td_dir(1).bin_size); 
                speed_data(end+1,:) = td_dir(i_trial).speed(idx_speed(1):idx_speed(2));
                
            end

            mean_speed_data(i_tgt,:) = mean(speed_data,1);
            
            tgt_dir_mean_fr(i_tgt) = mean(trial_fr);
            tgt_dir_std_fr(i_tgt) = std(trial_fr);
            
            figure(f_raster)
            subplot(3,3,subplot_map(i_tgt)); hold on
            plotRasterLIB(x_data,y_data,optsPlot,optsSave);
            
            figure(f_speed);
            ax_list(end+1) = subplot(3,3,subplot_map(i_tgt)); hold on
            plot(speed_x_data,mean_speed_data(i_tgt,:))
            xlim(window_plot)
        end
    catch
    end
% end

    linkaxes(ax_list,'xy');
%% fit tuning curve and plot

    cos_fit = fit(unique_tgt_dirs',tgt_dir_mean_fr','a*cos(x - b) + c','StartPoint',[15,0,15]);
    cos_x_data = linspace(-0.3927,2*pi-0.3927,1000);
    
    figure
    errorbar(unique_tgt_dirs*180/pi,tgt_dir_mean_fr,tgt_dir_std_fr,'k.','markersize',20)
    hold on
    plot(cos_x_data*180/pi,feval(cos_fit,cos_x_data),'k--')
    
    conf_vals = confint(cos_fit);
    plot(conf_vals(:,2)*180/pi + 360,[45,45],'k','linewidth',2)
    
    formatForLee(gcf)
    xlabel('Target direction (deg)');
    ylabel('Firing rate (Hz)');
    set(gca,'fontsize',14);
    
    xlim([-0.3927*180/pi,(2*pi-0.3927)*180/pi])
    
%% analyze movements during center hold caused by stim....
    window = [-20,100];
    td_stim = td_all(~isnan([td_all.idx_stimTime]) & [td_all.stimCode] >= 4);
    td_stim = trimTD(td_stim,{'idx_stimTime',window(1)},{'idx_stimTime',window(2)});
    x_data = [window(1):1:window(2)]*td_stim(1).bin_size;
    % get correction latency....
    % find max speed within 600 ms from stim onset, then go backwards to
    % get movement onset
    move_on_idx = nan(numel(td_stim),1);
    center_hold_speed = nan(numel(td_stim),1);
    
    color_idx = 1;
    figure();
    for i_trial = 1:numel(td_stim)
        [max_speed,max_speed_idx] = max(td_stim(i_trial).speed(td_stim(i_trial).idx_stimTime:end)); % relative to stim onset
        center_hold_speed(i_trial) = mean(td_stim(i_trial).speed(td_stim(i_trial).idx_stimTime-2:td_stim(i_trial).idx_stimTime+2));
        threshold = center_hold_speed(i_trial) + 0.1*(max_speed-center_hold_speed(i_trial));
        
        temp = find(td_stim(i_trial).speed(td_stim(i_trial).idx_stimTime:end) >= threshold,1,'first') - 1 + td_stim(i_trial).idx_stimTime - 1;
        if(~isempty(temp))
            move_on_idx(i_trial) = temp;
        end
    end
    
% make plots
    f_pos=figure(); hold on
    f_speed=figure(); hold on
    
    for i_trial = 1:numel(td_stim)
%         % pos plot
        if(i_trial < 15)
            figure(f_pos);
            xlabel('X-pos (cm)');
            ylabel('Y-pos (cm)');
            plot(td_stim(i_trial).pos(1:td_stim(i_trial).idx_stimTime,1),td_stim(i_trial).pos(1:td_stim(i_trial).idx_stimTime,2),...
                'color',getColorFromList(1,1))
             plot(td_stim(i_trial).pos(td_stim(i_trial).idx_stimTime:end,1),td_stim(i_trial).pos(td_stim(i_trial).idx_stimTime:end,2),...
                'color',getColorFromList(1,0))
        end
        %speed plot
        if(center_hold_speed(i_trial) < 5)
            figure(f_speed)

            subplot(3,1,1); hold on
            ylabel('X-vel (cm/s)');
            plot(x_data,td_stim(i_trial).vel(:,1),'color',getColorFromList(1,color_idx))
%             plot(x_data(move_on_idx(i_trial)),td_stim(i_trial).vel(move_on_idx(i_trial),1),'marker','.','markersize',20,'color',getColorFromList(1,2));
            
            subplot(3,1,2); hold on
            ylabel('Y-vel (cm/s)');
            plot(x_data,td_stim(i_trial).vel(:,2),'color',getColorFromList(1,color_idx))
%             plot(x_data(move_on_idx(i_trial)),td_stim(i_trial).vel(move_on_idx(i_trial),2),'marker','.','markersize',20,'color',getColorFromList(1,2));

            subplot(3,1,3); hold on
        hold on;
            ylabel('Speed (cm/s)');
            xlabel('Time after stimulation onset (s)');
            plot(x_data,td_stim(i_trial).speed,'color',getColorFromList(1,color_idx))
%             plot(x_data(move_on_idx(i_trial)),td_stim(i_trial).speed(move_on_idx(i_trial),1),'marker','.','markersize',20,'color',getColorFromList(1,2));
        end
    end


