%% set initial parameters

    input_data.folderpath = 'D:\Lab\Data\StimEvokedMove\Han_20200710_COBumpstim_100uA_singleChan_list1_250msTL\';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\6250-002085\SN 6250-002085.cmp';
    
    input_data.date = '20200610';
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
    stim_line_name = 'ainp16';
    
    file_name = dir('*nev*');
    
    params.event_list = {'goCueTime';'bumpTime';'bumpMagnitude';'stimTime';'stimCode';'isStimTrial';'isVisualTrial'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [0:1:255];
    params.trial_results = {'R','F','I','A'};
    
    move_onset_params.pre_move_thresh = 1000;
    move_onset_params.min_s = 3;
    move_onset_params.start_idx_offset = -10;
    move_onset_params.max_rt_offset = 400;
    
    td_all = [];
    for f = 1%:numel(file_name)
        cds = commonDataStructure();
        cds.file2cds(strcat(input_data.folderpath,file_name(f).name),input_data.array,input_data.monkey,input_data.ranBy,...
            input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
    
        
        % make trial data
        bad_trial_param.remove_nan_idx = true;
        bad_trial_param.nan_idx_names = {'idx_goCueTime','idx_endTime','idx_startTime'};
        td_temp = parseFileByTrial(cds,params);
        td_temp = getSpeed(td_temp);
        
        if(strcmpi(input_data.task,'taskRT'))
            td_temp = getGoCueTime(td_temp,cds);
            for i_trial = 1:numel(td_temp)
                td_temp(i_trial).idx_stimTime = td_temp(i_trial).idx_goCueTime;
            end
        end
        
        td_temp = removeBadTrials(td_temp,bad_trial_param);

        td_all = [td_all,td_temp];
    end
    
    if(td_all(1).bin_size < 0.05)
        % set it to 20ms
        td_all = binTD(td_all,ceil(0.01/td_all(1).bin_size));
    end
    
    
%% plot trajectory from before cue (stim or bump) to end of trial?

    % plot_mode = 1 for stim only, plot_mode = 2 for bump only, plot_mode =
    % 3 for stim and bump
    plot_mode = 1;
    cue_window = [-0.1,0.5]; % s
    cue_idx = ceil(cue_window/td_all(1).bin_size);
    unique_tgt_dirs = unique([td_all.target_direction]);
    
    switch plot_mode
        case 1
            td_plot = td_all(~isnan([td_all.stimCode]) & [td_all.bumpMagnitude] <=0);
            td_plot = trimTD(td_plot,{'idx_stimTime',cue_idx(1)},{'idx_stimTime',cue_idx(2)});
        case 2
            td_plot = td_all(isnan([td_all.stimCode]) & [td_all.bumpMagnitude] > 0);
            td_plot = trimTD(td_plot,{'idx_bumpTime',cue_idx(1)},{'idx_bumpTime',cue_idx(2)});
        case 3
            td_plot = td_all(~isnan([td_all.stimCode]) & [td_all.bumpMagnitude] > 0 & [td_all.target_direction] == unique_tgt_dirs(4));
            td_plot = trimTD(td_plot,{'idx_stimTime',cue_idx(1)},{'idx_stimTime',cue_idx(2)});
    end

    figure(); hold on;
    for i_trial = 1:numel(td_plot)
         plot(td_plot(i_trial).pos(:,1),td_plot(i_trial).pos(:,2),'k')
         
         if(plot_mode == 1)
             idx_plot = td_plot(i_trial).idx_stimTime;
             plot(td_plot(i_trial).pos(idx_plot,1),td_plot(i_trial).pos(idx_plot,2),'r.','markersize',12)

             idx_plot = td_plot(i_trial).idx_stimTime + ceil(0.1/td_plot(i_trial).bin_size);
             plot(td_plot(i_trial).pos(idx_plot,1),td_plot(i_trial).pos(idx_plot,2),'b.','markersize',12)

             idx_plot = td_plot(i_trial).idx_stimTime + ceil(0.20/td_plot(i_trial).bin_size);
             plot(td_plot(i_trial).pos(idx_plot,1),td_plot(i_trial).pos(idx_plot,2),'g.','markersize',12)
         elseif(plot_mode == 2)
             idx_plot = td_plot(i_trial).idx_bumpTime;
             plot(td_plot(i_trial).pos(idx_plot,1),td_plot(i_trial).pos(idx_plot,2),'r.','markersize',12)

             idx_plot = td_plot(i_trial).idx_bumpTime + ceil(0.1/td_plot(i_trial).bin_size);
             plot(td_plot(i_trial).pos(idx_plot,1),td_plot(i_trial).pos(idx_plot,2),'b.','markersize',12)

             idx_plot = td_plot(i_trial).idx_bumpTime + ceil(0.20/td_plot(i_trial).bin_size);
             plot(td_plot(i_trial).pos(idx_plot,1),td_plot(i_trial).pos(idx_plot,2),'g.','markersize',12)
         elseif(plot_mode == 3)
             if(~isnan(td_plot(i_trial).idx_stimTime) && ~isnan(td_plot(i_trial).idx_bumpTime))
                 idx_plot = td_plot(i_trial).idx_stimTime;
                 plot(td_plot(i_trial).pos(idx_plot,1),td_plot(i_trial).pos(idx_plot,2),'r.','markersize',12)

                 idx_plot = td_plot(i_trial).idx_stimTime + ceil(0.1/td_plot(i_trial).bin_size);
                 plot(td_plot(i_trial).pos(idx_plot,1),td_plot(i_trial).pos(idx_plot,2),'b.','markersize',12)

                 idx_plot = td_plot(i_trial).idx_bumpTime;
                 plot(td_plot(i_trial).pos(idx_plot,1),td_plot(i_trial).pos(idx_plot,2),'g.','markersize',12)
             end
         end
    end
    

%% plot speed on each trial for stim conditions, bump conditions and stim + bump conditions

    % plot_mode = 1 for stim only, plot_mode = 2 for bump only, plot_mode =
    % 3 for stim and bump
    plot_mode = 1;
    cue_window = [-0.1,0.3]; % s
    cue_idx = ceil(cue_window/td_all(1).bin_size);
    
    switch plot_mode
        case 1
            td_plot = td_all(~isnan([td_all.stimCode]) & [td_all.stimCode]==16 & ([td_all.bumpMagnitude] <=0 | isnan([td_all.idx_bumpTime])));
            td_plot = trimTD(td_plot,{'idx_stimTime',cue_idx(1)},{'idx_stimTime',cue_idx(2)});
        case 2
            td_plot = td_all(isnan([td_all.stimCode]) & [td_all.bumpMagnitude] > 0);
            td_plot = trimTD(td_plot,{'idx_bumpTime',cue_idx(1)},{'idx_bumpTime',cue_idx(2)});
        case 3
            td_plot = td_all(~isnan([td_all.stimCode]) & [td_all.bumpMagnitude] > 0);
            td_plot = trimTD(td_plot,{'idx_stimTime',cue_idx(1)},{'idx_stimTime',cue_idx(2)});
    end

    figure(); hold on;
    x_data = ((1:1:size(td_plot(1).speed))-1)*td_plot(1).bin_size + cue_window(1);
    for i_trial = 1:numel(td_plot)
         plot(x_data,td_plot(i_trial).speed,'k')
         
         
        if(plot_mode == 3)
             if(~isnan(td_plot(i_trial).idx_stimTime) && ~isnan(td_plot(i_trial).idx_bumpTime))
                 idx_plot = td_plot(i_trial).idx_bumpTime;
                 plot(x_data(idx_plot),td_plot(i_trial).speed(idx_plot),'g.','markersize',12)
             end
         end
    end





