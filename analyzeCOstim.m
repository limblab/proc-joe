%% set initial parameters

    input_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\CObump\Han_20200131_CObumpstim\';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\6250-002085\SN 6250-002085.cmp';
    
    input_data.date = '20200131';
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

    move_onset_params.pre_move_thresh = 1000;
    move_onset_params.min_s = 3;
    move_onset_params.start_idx_offset = -10;
    move_onset_params.max_rt_offset = 400;
    
    td_all = [];
    for f = 1:numel(file_name)
        cds = commonDataStructure();
        cds.file2cds(strcat(input_data.folderpath,file_name(1).name),input_data.array,input_data.monkey,input_data.ranBy,...
            input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
    
        % fix stim time based on analog input
        analog_idx = 1;
        for i = 1:numel(cds.analog)
            if(~isempty(find(strcmp(cds.analog{2}.Properties.VariableNames, stim_line_name))))
                analog_idx = i;
            end
        end
        
        stim_on = cds.analog{analog_idx}.t(find(diff(cds.analog{2}.(stim_line_name)-mean(cds.analog{2}.(stim_line_name))>3)>.5));
        
        stim_diff = nan(size(cds.trials,1),1);
        for tr = 1:size(cds.trials,1)
            if(~isnan(cds.trials.stimTime(tr)))
                stim_diff(tr,1) = stim_on(find(stim_on > cds.trials.stimTime(tr),1,'first')) - cds.trials.stimTime(tr);
            end
        end
        
        % make trial data
        bad_trial_param.remove_nan_idx = true;
        bad_trial_param.nan_idx_names = {'idx_goCueTime','idx_endTime','idx_startTime'};
        td_temp = parseFileByTrial(cds,params);
        td_temp = getSpeed(td_temp);
        
        td_temp = removeBadTrials(td_temp,bad_trial_param);

        td_all = [td_all,td_temp];
        for tr = 1:numel(td_all)
            if(~isnan(stim_diff(tr)))
                td_all(tr).idx_stimTime = td_all(tr).idx_stimTime + ceil(stim_diff(tr)*td_all(tr).bin_size);
            end
        end
    end
    
    if(td_all(1).bin_size < 0.01)
        % set it to 10ms
        td_all = binTD(td_all,ceil(0.01/td_all(1).bin_size));
    end


%% plot reaches for each direction from go cue to end
    offset = [0,0];
    subtract_start_point = 0;
    figure(); hold on;
    stim_codes = unique([td_all.stimCode]); stim_codes(isnan(stim_codes)) = [];
    
    for iCode = 1:numel(stim_codes)
        stim_lines = {};
        subplot(2,ceil(numel(stim_codes)/2),iCode); hold on;
        for tr = 1:numel(td_all)
            idx_goCue = td_all(tr).idx_goCueTime;
            idx_end = td_all(tr).idx_endTime;
            idx_stim = td_all(tr).idx_stimTime;

            color_plot = getColorFromList(1,1);
            if(~isnan(idx_stim))
                color_plot = getColorFromList(1,0);
            end

            if(isnan(td_all(tr).idx_bumpTime) && (isnan(td_all(tr).stimCode) || td_all(tr).stimCode == stim_codes(iCode)))
                x_data = td_all(tr).pos(idx_goCue+offset(1):idx_end+offset(2),1);
                y_data = td_all(tr).pos(idx_goCue+offset(1):idx_end+offset(2),2);
                if(subtract_start_point)
                    x_data = x_data - td_all(tr).pos(idx_goCue+offset(1),1);
                    y_data = y_data - td_all(tr).pos(idx_goCue+offset(2),2);
                end
                % plot line
                h=plot(x_data,y_data,...
                    '-','linewidth',1.5,'color',color_plot); hold on
                % plot stim time if applicable
                if(~isnan(idx_stim))
                    stim_lines{end+1} = h;
                    h=plot(x_data(idx_stim-(idx_goCue+offset(1))),y_data(idx_stim-(idx_goCue+offset(2))),...
                        '.','markersize',20,'color',getColorFromList(1,2));
                    stim_lines{end+1} = h;
                    h=plot(x_data(idx_stim-(idx_goCue+offset(1))+20),y_data(idx_stim-(idx_goCue+offset(2))+20),...
                        '.','markersize',20,'color',getColorFromList(1,4));
                    stim_lines{end+1} = h;
                end
                
            end
        end
        
        for i = 1:numel(stim_lines)
            uistack(stim_lines{i},'top');
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







    