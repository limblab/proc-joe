%% set initial parameters

    input_data.folderpath = 'D:\Lab\Data\CObumpmove\Han_20200616_stim\';
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
    
    params.event_list = {'goCueTime';'bumpTime';'stimTime';'stimCode'};
    params.extra_time = [1,2];
    params.include_ts = 0;
    params.exclude_units = [0:1:255];
    params.trial_results = {'R','A'};
    
    move_onset_params.pre_move_thresh = 1000;
    move_onset_params.min_s = 3;
    move_onset_params.start_idx_offset = -10;
    move_onset_params.max_rt_offset = 400;
    
    td_all = [];
    for f = 1:numel(file_name)
        cds = commonDataStructure();
        cds.file2cds(strcat(input_data.folderpath,file_name(f).name),input_data.array,input_data.monkey,input_data.ranBy,...
            input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
    
        
        % make trial data
        bad_trial_param.remove_nan_idx = true;
        bad_trial_param.nan_idx_names = {'idx_goCueTime','idx_endTime','idx_startTime'};
        td_temp = parseFileByTrial(cds,params);
        td_temp = getSpeed(td_temp);
        
        td_temp = removeBadTrials(td_temp,bad_trial_param);

        td_all = [td_all,td_temp];
    end
    
    if(td_all(1).bin_size < 0.01)
        % set it to 10ms
        td_all = binTD(td_all,ceil(0.01/td_all(1).bin_size));
    end


%% plot reaches for each direction from go cue to end
    offset = [0,75];
    subtract_start_point = 0;
    figure(); hold on;
    stim_codes = unique([td_all.stimCode]); stim_codes(isnan(stim_codes)) = [];
    max_plot = 100;
    
    for iCode = 1:numel(stim_codes)
        stim_lines = {};
        subplot(2,ceil(numel(stim_codes)/2),iCode); hold on;
        num_plot = 0;
        for tr = 1:numel(td_all)
            idx_goCue = td_all(tr).idx_goCueTime;
            idx_end = td_all(tr).idx_endTime;
            idx_stim = td_all(tr).idx_stimTime;

            color_plot = getColorFromList(1,1);
            if(~isnan(idx_stim))
                color_plot = getColorFromList(1,0);
            end

            if((isnan(td_all(tr).idx_bumpTime) && isnan(td_all(tr).stimCode) && num_plot < max_plot) || ... 
                    td_all(tr).stimCode == stim_codes(iCode))
                
                num_plot = num_plot + 1;
                x_data = td_all(tr).pos(idx_goCue+offset(1):idx_goCue+offset(2),1);
                y_data = td_all(tr).pos(idx_goCue+offset(1):idx_goCue+offset(2),2);
                if(subtract_start_point)
                    x_data = x_data - td_all(tr).pos(idx_goCue+offset(1),1);
                    y_data = y_data - td_all(tr).pos(idx_goCue+offset(1),2);
                end
                % plot line
                h=plot(x_data,y_data,...
                    '-','linewidth',1.5,'color',color_plot); hold on
                % plot stim time if applicable
                if(~isnan(idx_stim) && idx_stim > idx_goCue)
                    stim_lines{end+1} = h;
                    h=plot(x_data(idx_stim-(idx_goCue+offset(1))),y_data(idx_stim-(idx_goCue+offset(1))),...
                        '.','markersize',20,'color',getColorFromList(1,2));
                    stim_lines{end+1} = h;
                    h=plot(x_data(idx_stim-(idx_goCue+offset(1))+20),y_data(idx_stim-(idx_goCue+offset(1))+20),...
                        '.','markersize',20,'color',getColorFromList(1,4));
                    stim_lines{end+1} = h;
                end
                
            end
        end
        
        for i = 1:numel(stim_lines)
            uistack(stim_lines{i},'top');
        end
    end

%% plot average here
    stim_trial_mask = ~isnan([td_all.idx_stimTime]);
    no_bump_mask = isnan([td_all.idx_bumpTime]) | ([td_all.idx_bumpTime] < [td_all.idx_goCueTime]);
    
    td_stim = td_all(stim_trial_mask == 1);
    td_no_stim = td_all(stim_trial_mask == 0 & no_bump_mask == 0);
    
    td_stim = trimTD(td_stim,{'idx_stimTime',-20},{'idx_stimTime',45});
    td_no_stim = trimTD(td_no_stim,{'idx_goCueTime',0},{'idx_goCueTime',55});

    td_no_stim_avg = trialAverage(td_no_stim,'target_direction');
    
    figure(); hold on;
    for i_code = 1:numel(stim_codes)
        try
            td_stim_code = td_stim([td_stim.stimCode] == stim_codes(i_code));
            td_stim_avg = trialAverage(td_stim_code,'target_direction');
            subplot(2,4,i_code)
            for i_trial = 1:numel(td_no_stim_avg)
                idx_stim = td_stim_avg(i_trial).idx_stimTime;
                idx_stim_offset = idx_stim + floor(0.2/td_stim_avg(i_trial).bin_size);
                plot(td_stim_avg(i_trial).pos(:,1),td_stim_avg(i_trial).pos(:,2),'color',getColorFromList(1,0)); hold on
                plot(td_no_stim_avg(i_trial).pos(:,1),td_no_stim_avg(i_trial).pos(:,2),...
                    'color',getColorFromList(1,1)); hold on
                plot(td_stim_avg(i_trial).pos(idx_stim,1),td_stim_avg(i_trial).pos(idx_stim,2),...
                    '.','markersize',20,'color',getColorFromList(1,2)); hold on
                plot(td_stim_avg(i_trial).pos(idx_stim_offset,1),td_stim_avg(i_trial).pos(idx_stim_offset,2),...
                    '.','markersize',20,'color',getColorFromList(1,4)); hold on

            end
        end
    end
    
%% plot PSTH for each reach direction and during CO for stim case for each unit
    array_name = input_data.array(6:end);
    spike_field = [array_name,'_spikes'];
    unit_guide_field = [array_name,'_unit_guide'];
    
    % separate trials into stim center hold, stim move, and no stim
    stim_trial_mask = ~isnan([td_all.idx_stimTime]);
    td_stim = td_all(stim_trial_mask == 1);
    td_no_stim = td_all(stim_trial_mask == 0);
    
    stim_hold_mask = [td_stim.idx_stimTime] < [td_stim.idx_goCueTime];
    td_stim_hold = td_stim(stim_hold_mask == 1);
    td_stim_move = td_stim(stim_hold_mask == 0);
    
    % trim TD stim around stimulation
    x_data_stim = [-7:12];
    x_data_hold = [-27:-8];
    x_data_move = [-3:16];
    
    td_stim_hold = trimTD(td_stim_hold,{'idx_stimTime',x_data_stim(1)},{'idx_stimTime',x_data_stim(end)});
    td_stim_move = trimTD(td_stim_move,{'idx_stimTime',x_data_stim(1)},{'idx_stimTime',x_data_stim(end)});
    % trim TD no stim around center hold and movement
    td_no_hold = trimTD(td_no_stim,{'idx_goCueTime',x_data_hold(1)},{'idx_goCueTime',x_data_hold(end)});
    td_no_move = trimTD(td_no_stim,{'idx_goCueTime',x_data_move(1)},{'idx_goCueTime',x_data_move(end)});
    
    % make PSTH for each reach direction and neuron
    tgt_dirs = unique([td_stim_move.target_direction]);
    subplot_idx = [6,2,8,4];
    for nn = 1:size(td_all(1).(spike_field),2)
        figure();
        for tgt_idx = 1:numel(tgt_dirs)
            tgt_dir_mask = [td_stim_move.target_direction] == tgt_dirs(tgt_idx);
            td_stim_tgt_dir = td_stim_move(tgt_dir_mask == 1);
            
            sig = squeeze(getSigByTrial(td_stim_tgt_dir,{spike_field,nn}));
            % sig is an time point x trial matrix. Average over trials
            sig_mean = mean(sig,2);
            
            % plot sig_mean on correct subplot
            subplot(3,3,subplot_idx(tgt_idx))
            
            plot(x_data_stim,sig_mean,'r');
            hold on
            % do same for no stim trials
            tgt_dir_mask = [td_no_move.target_direction] == tgt_dirs(tgt_idx);
            td_no_tgt_dir = td_no_move(tgt_dir_mask == 1);
            
            sig = squeeze(getSigByTrial(td_no_tgt_dir,{spike_field,nn}));
            % sig is an time point x trial matrix. Average over trials
            sig_mean = mean(sig,2);
            plot(x_data_stim,sig_mean,'k');
            
            xlim([x_data_stim(1),x_data_stim(end)]);
        end
        
        % plot center hold data
        sig = squeeze(getSigByTrial(td_stim_hold,{spike_field,nn}));
        % sig is an time point x trial matrix. Average over trials
        sig_mean = mean(sig,2);

        % plot sig_mean on correct subplot
        subplot(3,3,5)

        plot(x_data_stim,sig_mean,'r');
        hold on;
        % do same for no stim data
        sig = squeeze(getSigByTrial(td_no_hold,{spike_field,nn}));
        % sig is an time point x trial matrix. Average over trials
        sig_mean = mean(sig,2);
        
        plot(x_data_stim,sig_mean,'k');
        xlim([x_data_stim(1),x_data_stim(end)]);
    end







    