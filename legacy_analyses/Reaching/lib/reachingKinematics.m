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
            f_hand_handle.Name = [td_list{i_td}.monkey, '_', td_list{i_td}.date, '_', task_list{i_td}, '_handVsHandlePos'];
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
%   Compute workspace volume for each task, and compute x,y,z length
%   Plot elbow position against each other in a 3D plot
    
    f_workspace(1) = figure();
    f_workspace(2) = figure();
    f_corr = figure();
    task_3d_idx = find(strcmpi(task_list,'RT3D')|strcmpi(task_list,'freeReach'),1,'first');
    task_2d_idx = find(strcmpi(task_list,'RT'),1,'first');
    
    markernames = {'hand2','elbow1'};
    perc_bound = [5,95];
    
    workspace_vol = nan(1,2);
    workspace_area = nan(1,2);
    workspace_min = nan(3,2); % x,y,z; 2D and 3D task
    workspace_max = nan(3,2);
    
    corr_V_RT2D = [];
    corr_S_RT2D = [];
    corr_V_RT3D = [];
    corr_S_RT3D = [];
    
    kin_corr = []; % corr matrix for each task (will be a 3D array)
    for i_task = 1:2 %loop for both RT2D and RT3D tasks
        switch i_task
            case 1
                task_idx = task_2d_idx;
            case 2
                task_idx = task_3d_idx;
        end
        
        %start of plotting hand marker distribution, and correlation
        %between markers
        for i_marker = 1:numel(markernames) %loop for each marker
            % plot hand position for each task
            figure(f_workspace(i_marker));
            f_workspace(i_marker).Name = [td_list{task_idx}.monkey, '_', td_list{task_idx}.date, '_', task_list{task_idx}, '_',markernames{i_marker},'Pos'];
            
            markername = markernames{i_marker};
            title(markername)
            dlc_idx = [find((strcmpi(td_list{task_idx}.dlc_pos_names,[markername,'_x']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,[markername,'_y']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,[markername,'_z'])))];
            
            temp_start = 1;
            temp_end = length(td_list{task_idx}.dlc_pos);
            
            plot3(td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(1)),td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(2)),td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(3)),...
                '.','linestyle','none','color',getColorFromList(1,i_task-1)); hold on;

            grid on
            xlabel('x-pos (cm)'); ylabel('y-pos (cm)'); zlabel('z-pos');
            l=legend('RT2D','RT3D');
            set(l,'box','off','location','best');
            formatForLee(gcf);
            set(gca,'fontsize',14);
            
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
        %end of plotting hand marker distribution, and correlation
        %between markers
        
        
        
        % for 2D and 3D task, compute correlation between different arm
        % markers
        corr_markernames = {'elbow1','hand2'};
        dlc_idx = [];
        for i_marker = 1:numel(corr_markernames)
            dlc_idx = [dlc_idx,find((strcmpi(td_list{task_idx}.dlc_pos_names,[corr_markernames{i_marker},'_x']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,[corr_markernames{i_marker},'_y']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,[corr_markernames{i_marker},'_z'])))];
        end
        
        % get velocity for each marker
        dlc_data_vel = td_list{task_idx}.dlc_vel(:,dlc_idx);
        % remove nan's
        dlc_data_vel = dlc_data_vel(~any(isnan(dlc_data_vel),2),:);
        %Calculate speed
        dlc_data_speed = [sqrt(sum(dlc_data_vel(:,1:2).^2,2)),sqrt(sum(dlc_data_vel(:,3:4).^2,2))];        
        
        corr_vel = corr(dlc_data_vel);
        corr_speed = corr(dlc_data_speed);
        
        %Save the correlation data for a scatterplot
        if task_idx == task_3d_idx
            corr_V_RT3D = corr_vel;
            corr_S_RT3D = corr_speed;
        else
            corr_V_RT2D = corr_vel;
            corr_S_RT2D = corr_speed;
        end
    end
  
    %Make a scatter plot with the correlation data
    figure(f_corr); hold on;
    f_corr.Name = [td_list{task_idx}.monkey, '_', td_list{task_idx}.date, '_', task_list{task_idx}, '_corr2Dvs3D'];
    % plot hand x (1) vs elbow x (4)
    plot(corr_V_RT2D(1,4),corr_V_RT3D(1,4),'r.','markersize',40);
    % plot hand y (2) vs elbow y (5)
    plot(corr_V_RT2D(2,5),corr_V_RT3D(2,5),'c.','markersize',40);
    % plot hand z (3) vs elbow z(6)
    plot(corr_V_RT2D(3,6), corr_V_RT3D(3,6),'m.','markersize',40);
    % plot hand speed (1) vs elbow speed (2)
    plot(corr_S_RT2D(1,2),corr_S_RT3D(1,2),'k.','markersize',40);
   
    kin_corr.corr_3D = [corr_V_RT3D(1,4), corr_V_RT3D(2,5), corr_V_RT3D(3,6), corr_S_RT3D(1,2)];
    kin_corr.corr_2D = [corr_V_RT2D(1,4), corr_V_RT2D(2,5), corr_V_RT2D(3,6), corr_S_RT2D(1,2)];
    
    
    %Set the lines on the graph
    h(1) = plot([0 0],[-1000 1000],'k-','linewidth',0.5); % plot y-axis
    h(2) = plot([-1000 1000],[0 0],'k-','linewidth',0.5); % plot x-axis
    h(3) = plot([-1000 1000],[-1000 1000],'k--','linewidth',0.5); % plot unity line
    set( get( get( h(1), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    set( get( get( h(2), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    set( get( get( h(3), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    lgd = legend('Elbow V_X - Hand V_X',...
    'Elbow V_Y - Hand V_Y',...
    'Elbow V_Z - Hand V_Z',...
    'Elbow spd - Hand spd',...
    '','','');
    lgd.Location = 'best';
    xlim([-0.3,1]);
    ylim([-0.3,1]);
    xlabel('RT2D Correlation');
    ylabel('RT3D Correlation');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    
    % plot correlation between joint angles for each task as a heatmap,
    % include difference between tasks
    
    f_joint = [];
    corr_joint = {};
    max_corr = -1;
    for i_td = 1:numel(td_list)
        if(isfield(td_list{i_td},'opensim_names') && ~isempty(td_list{i_td}.opensim_names))
            markername = 'vel';
            [vel_idx] = find(~cellfun(@isempty,strfind(td_list{i_td}.opensim_names,strcat('_',markername))));

            joint_vel = td_list{i_td}.opensim(:,vel_idx);
            joint_vel_names = erase(td_list{i_td}.opensim_names(vel_idx),'_vel'); % remove trailing identifier
            joint_vel_names = strrep(joint_vel_names,'_',' ');
            
            corr_joint{i_td} = corr(joint_vel);
            for i_joint = 1:numel(joint_vel_names) % overwrite diagonal vals
                corr_joint{i_td}(i_joint,i_joint) = nan;
            end
            max_corr = max(max_corr,max(max(corr_joint{i_td},[],'omitnan'),[],'omitnan'));
        end
    end
    diff_corr = abs(corr_joint{task_3d_idx})-abs(corr_joint{task_2d_idx});
    max_corr = max(max_corr,abs(max(max(diff_corr,[],'omitnan'),[],'omitnan')));
    
    % make figs
    im_alpha = tril(ones(numel(joint_vel_names)),-1); % set diagonal and upper diagonal as not visible
    for i_td = 1:numel(td_list)
        f_joint = figure();
        f_joint.Name = [td_list{i_td}.monkey, '_', td_list{i_td}.date, '_',td_list{i_td}.task, '_jointAngleCorrelation'];
        f_joint.Position = [691 124 1099 839];
        
        imagesc(abs(corr_joint{i_td}),'AlphaData',im_alpha);
        caxis([0,ceil(max_corr*10)/10]); % round up to nearest tenth
        colormap(colorcet('L5'));
        b=colorbar;
        set(gca,'XTick',1:1:numel(joint_vel_names),'XTickLabel',joint_vel_names);
        set(gca,'YTick',1:1:numel(joint_vel_names),'YTickLabel',joint_vel_names);
        
    end
    % plot difference between correlations
    f_joint = figure();
    f_joint.Name = [td_list{i_td}.monkey, '_', td_list{i_td}.date, '_jointAngleCorrelationDifference'];
    f_joint.Position = [691 124 1099 839];
    
    imagesc(diff_corr,'AlphaData',im_alpha);
    caxis([-ceil(max_corr*10)/10,ceil(max_corr*10)/10]); % round up to nearest tenth
    colormap(colorcet('D4'));
    b=colorbar;
    set(gca,'XTick',1:1:numel(joint_vel_names),'XTickLabel',joint_vel_names);
    set(gca,'YTick',1:1:numel(joint_vel_names),'YTickLabel',joint_vel_names);
    
    
    % plot speed histogram for hand and elbow markers, comparing both tasks
    markernames = {'hand2','elbow1'};
    spd_data = cell(numel(td_list),1); % speed data for each trial data and marker
    vel_data = cell(numel(td_list),1); % vel data for each trial data and marker
    max_spd = 0;
    max_vel = 0;
    for i_td = 1:numel(td_list)
        for i_marker = 1:numel(markernames)
            dlc_idx = [find((strcmpi(td_list{i_td}.dlc_pos_names,[markernames{i_marker},'_x']))),...
                find((strcmpi(td_list{i_td}.dlc_pos_names,[markernames{i_marker},'_y']))),...
                find((strcmpi(td_list{i_td}.dlc_pos_names,[markernames{i_marker},'_z'])))];
            
            marker_vel = td_list{i_td}.dlc_vel(:,dlc_idx);
            marker_spd = sqrt(sum(marker_vel.^2,2));
            
            spd_data{i_td} = [spd_data{i_td}, marker_spd]; % append for each marker
            vel_data{i_td} = [vel_data{i_td}, marker_vel]; % append for each marker, x then y then z
            max_spd = max(max_spd, max(marker_spd));
            max_vel = max(max_vel, max(max(abs(marker_vel))));
        end
    end
    
    f_speed = figure();
    f_speed.Name = [td_list{i_td}.monkey, '_', td_list{i_td}.date, '_handElbowSpeedHist'];
    ax_list = [];
    bin_edges = linspace(0,max_spd,100);
    for i_marker = 1:numel(markernames)
        ax_list(end+1) = subplot(numel(markernames),1,i_marker); hold on
        for i_td = 1:numel(spd_data)
            switch i_td
                case 1
                    task_idx = task_2d_idx;
                case 2
                    task_idx = task_3d_idx;
            end
            
            [spd_counts] = histcounts(spd_data{task_idx}(:,i_marker),bin_edges,'normalization','probability');
            histogram('BinEdges',bin_edges,'BinCounts',spd_counts,'FaceColor','none','edgeColor',getColorFromList(1,i_td-1),'linewidth',2);
        end
        xlabel([markernames{i_marker},' speed (cm/s)']);
        ylabel('Proportion of data');
        formatForLee(gcf);
        set(gca,'fontsize',14);
    end
    linkaxes(ax_list,'xy');
    l=legend(task_list); set(l,'box','off');
    
    
    f_vel_dist = figure();
    f_vel_dist.Name = [td_list{i_td}.monkey, '_', td_list{i_td}.date, '_handElbowVelHist'];
    ax_list = [];
    bin_edges = linspace(-max_vel,max_vel,200);
    vel_names = {'x','y','z'};
    for i_marker = 1:numel(markernames)
        for i_dir = 1:3 % x, y, z
            ax_list(end+1) = subplot(3,numel(markernames),i_marker + 2*(i_dir-1)); hold on
            for i_td = 1:numel(vel_data)
                [vel_counts] = histcounts(vel_data{i_td}(:,i_dir + 3*(i_marker-1)),bin_edges,'normalization','probability');
                histogram('BinEdges',bin_edges,'BinCounts',vel_counts,'FaceColor','none','edgeColor',getColorFromList(1,i_td-1),'linewidth',2);
            end
            xlabel([markernames{i_marker},' Vel ' vel_names{i_dir} ' (cm/s)']);
            ylabel('Proportion of data');
            formatForLee(gcf);
            set(gca,'fontsize',14);
        end
    end
    linkaxes(ax_list,'xy');
    l=legend(task_list); set(l,'box','off');
    xlim([-max_vel,max_vel])
    % package outputs into nice tables so the data is easily worked with
    output_data.kin_corr = kin_corr;
    output_data.workspace_min = workspace_min;
    output_data.workspace_max = workspace_max;
    output_data.workspace_extent = workspace_max-workspace_min;
    output_data.task_list = {task_list{task_3d_idx},task_list{task_2d_idx}};
    output_data.workspace_vol = workspace_vol;
    output_data.workspace_area = workspace_area;

end

function [speed] = calculateSpeed(posArr, binSize)
    %velArr = diff(posArr/binSize);
    x2 = diff(posArr(:,1)).^2;
    y2 = diff(posArr(:,2)).^2;
    z2 = diff(posArr(:,3)).^2;
    speed = sqrt(x2 + y2 + z2)./binSize./1000; %mm/0.05second bin to> m/s
    %speed = sqrt(velArr(:,1).^2 + velArr(:,2).^2 + velArr(:,3).^2)/binSize;
    speed(end+1,:) = speed(end,:);
end