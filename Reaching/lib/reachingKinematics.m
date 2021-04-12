function [output_data] = reachingKinematics(td_list,task_list,input_data)
% for RT2D and RT3D tasks (in task list), this function gets and plots some relevant
% kinematic data:
% for 2D task, plot hand position against handle position
% for 2D and 3D task, plot hand position against each other in a 3D plot.
%   Compute workspace volume for each task and compares, also computes x,
%   y and z length
%   Plot elbow position against each other in a 3D plot
% for 2D and 3D task, compute correlation between arm kinematics
    output_data = [];


% for 2D task, plot hand position against handle position
    f_hand_handle = figure();
    num_2D_task = sum(strcmpi(task_list,'RT'));
    sub_counter = 1;
    for i_td = 1:numel(td_list)
        if(strcmpi(task_list{i_td},'RT')==1) % only plot for 2D task
            
            markername = 'hand2';
            dlc_idx = [find((strcmpi(td_list{1}.dlc_pos_names,[markername,'_x']))),find((strcmpi(td_list{1}.dlc_pos_names,[markername,'_y'])))];

            mean_td_pos = mean(td_list{i_td}.pos);
            mean_dlc_pos = mean(td_list{i_td}.dlc_pos);
            
            x_data = ((1:1:size(td_list{i_td}.pos,1))-1)*td_list{i_td}.bin_size;
            ax1=subplot(2,num_2D_task,sub_counter);
            plot(x_data,td_list{i_td}.pos(:,1)-mean_td_pos(1),'color',getColorFromList(1,0),'linestyle','-')
            hold on;
            plot(x_data,(td_list{i_td}.dlc_pos(:,dlc_idx(1))-mean_dlc_pos(dlc_idx(1))),'color',getColorFromList(1,1),'linestyle','-');
            
            ylabel('x-pos (cm)');
            l=legend('handle','dlc');
            set(l,'box','off','location','best');
            title('hand position against handle position');
            formatForLee(gca);
            set(gca,'fontsize',14);
            
            ax2=subplot(2,num_2D_task,sub_counter+num_2D_task);
            plot(x_data,td_list{i_td}.pos(:,2)-mean_td_pos(2),'color',getColorFromList(1,0),'linestyle','-')
            hold on
            plot(x_data,td_list{i_td}.dlc_pos(:,dlc_idx(2))-mean_dlc_pos(dlc_idx(2)),'color',getColorFromList(1,1),'linestyle','-');
            
            xlabel('Experiment time (s)');
            ylabel('y-pos (s)');
            formatForLee(gca);
            set(gca,'fontsize',14);
            
            sub_counter = sub_counter + 1;
            linkaxes([ax1,ax2],'x');
        end
    end

% for 2D and 3D task, plot hand position against each other in a 3D plot.
%   Compute workspace volume for each task and compares, also computes x,
%   y and z length
%   Plot elbow position against each other in a 3D plot
    
    f_workspace(1) = figure();
    f_workspace(2) = figure();
    f_corr = figure();
    task_3d_idx = find(strcmpi(task_list,'RT3D'),1,'first');
    task_2d_idx = find(strcmpi(task_list,'RT'),1,'first');
    markernames = {'hand2','elbow1'};
    perc_bound = [5,95];
    
    workspace_vol = nan(1,2);
    workspace_area = nan(1,2);
    workspace_min = nan(3,2); % x,y,z; 2D and 3D task
    workspace_max = nan(3,2);
    
    kin_corr = []; % corr matrix for each task (will be a 3D array)
    for i_task = 1:2
        switch i_task
            case 1
                task_idx = task_3d_idx;
                title_str = 'RT3D';
            case 2
                task_idx = task_2d_idx;
                title_str = 'RT2D';
        end
        
        for i_marker = 1:numel(markernames)
            % plot hand position for each task
            figure(f_workspace(i_marker));
            markername = markernames{i_marker};
            title(markername)
            dlc_idx = [find((strcmpi(td_list{task_idx}.dlc_pos_names,[markername,'_x']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,[markername,'_y']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,[markername,'_z'])))];
            
            temp_start = 5000;
            temp_end = 5200;
            
            %plot3(td_list{task_idx}.dlc_pos(:,dlc_idx(1)),td_list{task_idx}.dlc_pos(:,dlc_idx(2)),td_list{task_idx}.dlc_pos(:,dlc_idx(3)),...
            %    '.','linestyle','none','color',getColorFromList(1,i_task-1)); hold on;
            %plot3(td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(1)),td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(2)),td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(3)),...
            %    '.','linestyle','none','color',getColorFromList(1,i_task-1)); hold on;
            plot3(td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(1)),td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(2)),td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(3)),...
                'color',getColorFromList(1,i_task-1)); hold on;
            grid on
            xlabel('x-pos (cm)'); ylabel('y-pos (cm)'); zlabel('z-pos');
            l=legend('RT3D','RT2D');
            set(l,'box','off','location','best');
            
            % compute volume, and extent along each axis
            if(~isempty(strfind(markername,'hand')))
                % for extent, only look at data within some percentile
                bounds = prctile(td_list{task_idx}.dlc_pos(:,dlc_idx(:)),perc_bound);
                workspace_min(:,i_task) = min(bounds)';
                workspace_max(:,i_task) = max(bounds)';
                % then use only that data to compute the convhull
                data = td_list{task_idx}.dlc_pos(:,dlc_idx(:));
                data = data(all(data < workspace_max(:,i_task)' & data > workspace_min(:,i_task)',2),:);
                
                [~,workspace_vol(i_task)] = convhull(data(:,1),data(:,2),data(:,3));
                [~,workspace_area(i_task)] = convhull(data(:,1),data(:,2));
            end
        end

        
        % for 2D and 3D task, compute correlation between different arm
        % markers
        figure(f_corr)
        corr_markernames = {'shoulder','elbow1','hand2'};
        dlc_idx = [];
        for i_marker = 1:numel(corr_markernames)
            dlc_idx = [dlc_idx,find((strcmpi(td_list{task_idx}.dlc_pos_names,[corr_markernames{i_marker},'_x']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,[corr_markernames{i_marker},'_y']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,[corr_markernames{i_marker},'_z'])))];
        end
        dlc_data = td_list{task_idx}.dlc_vel(:,dlc_idx);
        dlc_data = dlc_data(~any(isnan(dlc_data),2),:);
                
        corr_data = corr(dlc_data);
        subplot(1,3,i_task)
        imagesc(corr_data,[-1,1]);
        set(gca, 'xTick', [1:9]);
        b=colorbar();
        kin_corr(:,:,i_task) = corr_data;
        title(title_str);
    end
     
    % plot correlation difference for 2D and 3D tasks, format figure
    figure(f_corr)
    subplot(1,3,3)
    imagesc(kin_corr(:,:,1) - kin_corr(:,:,2),[-1,1]);
    set(gca, 'xTick', [1:9]);
    b=colorbar();
    title('3D - 2D');
    
    % format workspace figures
    for i = 1:numel(f_workspace)
        figure(f_workspace(i))
        l=legend('3D task','2D task');
        set(l,'box','off','location','best');
        formatForLee(gcf);
        set(gca,'fontsize',14);
    end
    



    % package outputs into nice tables so the data is easily worked with
    output_data.kin_corr = kin_corr;
    output_data.workspace_min = workspace_min;
    output_data.workspace_max = workspace_max;
    output_data.workspace_extent = workspace_max-workspace_min;
    output_data.task_list = {task_list{task_3d_idx},task_list{task_2d_idx}};
    output_data.workspace_vol = workspace_vol;
    output_data.workspace_area = workspace_area;

end