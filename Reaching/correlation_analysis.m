%% load in TD_list

%% clean up data
% resample trial data to appropriate bin size
    bin_size = 0.05; % s
    for i_td = 1:numel(td_list)
%         if(td_list{i_td}.bin_size <= bin_size)
%             td_list{i_td} = binTD(td_list{i_td},bin_size/td_list{i_td}.bin_size);
%         else
%             warning('td bin size is larger than desired bin size');
%         end


        if(isfield(td_list{i_td},'dlc_pos'))
            % set origin as shoulder position at t=0
            td_list{i_td} = setOriginAsShoulder(td_list{i_td},1);
            % get marker velocity
            td_list{i_td} = getDifferential(td_list{i_td},struct('signals','dlc_pos','alias','dlc_vel'));
            % remove time points where dlc tracking is bad
            dlc_idx = [find((strcmpi(td_list{i_td}.dlc_pos_names,['hand2','_x']))),...
                        find((strcmpi(td_list{i_td}.dlc_pos_names,['hand2','_y']))),...
                        find((strcmpi(td_list{i_td}.dlc_pos_names,['elbow1','_x']))),...
                        find((strcmpi(td_list{i_td}.dlc_pos_names,['elbow1','_y'])))]; % not clear to me why we do x and y since 1 should be sufficient
            bad_points = any(isnan(td_list{i_td}.dlc_pos(:,dlc_idx)),2) | any(isnan(td_list{i_td}.dlc_vel(:,dlc_idx)),2);
            td_names = fieldnames(td_list{i_td});
            for i_name = 1:numel(td_names)
                if(size(bad_points,1) == size(td_list{i_td}.(td_names{i_name}),1))
                    td_list{i_td}.(td_names{i_name})(bad_points,:) = [];
                end
            end
            fprintf('Removed %d percent of trials because of missing markers\n',sum(bad_points)/numel(bad_points)*100)
            
            % fix x and y axes for dlc_pos dlc_vel and hand and elbow only
            % elbow
%             temp_dlc_x_pos = td_list{i_td}.dlc_pos(:,dlc_idx(3));
%             temp_dlc_x_vel = td_list{i_td}.dlc_vel(:,dlc_idx(3));
%             td_list{i_td}.dlc_pos(:,dlc_idx(3)) = -td_list{i_td}.dlc_pos(:,dlc_idx(4));
%             td_list{i_td}.dlc_vel(:,dlc_idx(3)) = -td_list{i_td}.dlc_vel(:,dlc_idx(4));
%             td_list{i_td}.dlc_pos(:,dlc_idx(4)) = temp_dlc_x_pos;
%             td_list{i_td}.dlc_vel(:,dlc_idx(4)) = temp_dlc_x_vel;
%             
%             % hand
%             temp_dlc_x_pos = td_list{i_td}.dlc_pos(:,dlc_idx(1));
%             temp_dlc_x_vel = td_list{i_td}.dlc_vel(:,dlc_idx(1));
%             td_list{i_td}.dlc_pos(:,dlc_idx(1)) = -td_list{i_td}.dlc_pos(:,dlc_idx(2));
%             td_list{i_td}.dlc_vel(:,dlc_idx(1)) = -td_list{i_td}.dlc_vel(:,dlc_idx(2));
%             td_list{i_td}.dlc_pos(:,dlc_idx(2)) = temp_dlc_x_pos;
%             td_list{i_td}.dlc_vel(:,dlc_idx(2)) = temp_dlc_x_vel;
%             
        end
        
        
        
        % remove trials where monkey's arm position is out of the
        % "workspace"
        td_list{i_td} = getExperimentPhase(td_list{i_td},task_list{i_td});
    end

%% compute cross-correlation in sliding window for each trial data

    window_size = 500;
    window_step = 500;
    figure(); hold on;
    corr_list = {};
    for i_td = 1:numel(td_list)
        td_temp = td_list{i_td};
        dlc_idx = [find((strcmpi(td_temp.dlc_pos_names,['hand2','_x']))),...
                        find((strcmpi(td_temp.dlc_pos_names,['hand2','_y']))),...
                        find((strcmpi(td_temp.dlc_pos_names,['elbow1','_x']))),...
                        find((strcmpi(td_temp.dlc_pos_names,['elbow1','_y'])))];
                    
        window_starts = 1:window_step:(size(td_temp.dlc_pos,1)-window_size+1);
        
        corrs = zeros(size(window_starts,1),3);
        
        for i_win = 1:numel(window_starts)
            data = td_temp.dlc_vel(window_starts(i_win):window_starts(i_win)+window_size-1, dlc_idx);
            [r,lags] = xcorr(data(:,2),data(:,4),50,'biased');
            plot(lags,r);
            % get elbow-x hand-x corr
            corrs(i_win,1) = corr(data(:,3),data(:,1));
            % get elbow-y hand-y corr
            corrs(i_win,2) = corr(data(:,4),data(:,2));

        end
        
        corr_list{i_td} = corrs;
        
    end


    %%
    task_2d = find(cellfun(@isempty,strfind(task_list,'3D')),1,'first');
    task_3d = find(~cellfun(@isempty,strfind(task_list,'3D')),1,'first');
    figure(); hold on;
    mean_corr = cellfun(@abs,cellfun(@mean,corr_list,'UniformOutput',false),'UniformOutput',false);
    std_corr = cellfun(@std,corr_list,'UniformOutput',false);
    
        % plot elbow-x hand-x corr
        idx = 1;
        errorbar(mean_corr{task_2d}(:,idx),mean_corr{task_3d}(:,idx),...
            std_corr{task_3d}(:,idx), std_corr{task_3d}(:,idx),...
            std_corr{task_2d}(:,idx), std_corr{task_2d}(:,idx),'.','markersize',20)
        
        % plot elbow-y hand-y corr
        idx = 2;
        errorbar(mean_corr{task_2d}(:,idx),mean_corr{task_3d}(:,idx),...
            std_corr{task_3d}(:,idx), std_corr{task_3d}(:,idx),...
            std_corr{task_2d}(:,idx), std_corr{task_2d}(:,idx),'.','markersize',20)
        
        % plot elbow speed hand speed corr


    plot([0,1],[0,1],'k--')
    