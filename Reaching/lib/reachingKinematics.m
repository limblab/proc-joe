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
            %plot(x_data,(td_list{i_td}.dlc_pos(:,dlc_idx(2))-mean_dlc_pos(dlc_idx(2))),'color',getColorFromList(1,1),'linestyle','-');
            
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
            %plot(x_data,td_list{i_td}.dlc_pos(:,dlc_idx(1))-mean_dlc_pos(dlc_idx(1)),'color',getColorFromList(1,1),'linestyle','-');
            
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
    task_3d_idx = find(strcmpi(task_list,'RT3D'),1,'first')
    task_2d_idx = find(strcmpi(task_list,'RT'),1,'first')
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
                task_idx = task_3d_idx;
                title_str = 'RT3D';
            case 2
                task_idx = task_2d_idx;
                title_str = 'RT2D';
        end
        
        %Plot a few consecutive frames on the whole arm for either RT2D and RT3D task
        %shoulder:  1,2,3
        %elbow2:    7,8,9
        %wrist2:    13,14,15
        %hand3:     22,23,24
        %task:      td_list{i_task}
        %pos data:  td_list{i_task}.dlc_pos(wa_f_st:wa_f_end , 1:3)
        %kinematics data: 
        wa_f_st = 5000;                 %whole arm frame start
        wa_f_len = 30;                  %whole arm frame length (frames)
        wa_f_end = wa_f_st + wa_f_len; %whole arm frame end
        marker_size = 75;
        marker_width = 3;
        line_width = 0.5;
        
        shoulder_col = 1:3;
        elbow2_col = 7:9;
        wrist2_col = 13:15;
        hand3_col = 22:24;
        
        one_in_n_frames = 1;
        
        shoulder_pos = td_list{i_task}.dlc_pos(wa_f_st:one_in_n_frames:wa_f_end , shoulder_col);
        elbow2_pos = td_list{i_task}.dlc_pos(wa_f_st:one_in_n_frames:wa_f_end , elbow2_col);
        wrist2_pos = td_list{i_task}.dlc_pos(wa_f_st:one_in_n_frames:wa_f_end , wrist2_col);
        hand3_pos = td_list{i_task}.dlc_pos(wa_f_st:one_in_n_frames:wa_f_end , hand3_col);

        shoulder_spd = calculateSpeed(shoulder_pos,td_list{i_task}.bin_size);
        elbow2_spd = calculateSpeed(elbow2_pos,td_list{i_task}.bin_size);
        wrist2_spd = calculateSpeed(wrist2_pos,td_list{i_task}.bin_size);
        hand3_spd = calculateSpeed(hand3_pos,td_list{i_task}.bin_size);

        %since the hand speed is roughly less than 2m/s, I'm just gonna set
        %the color map to 200 colors, and plot the colors accordingly
        shoulder_clr = round(shoulder_spd.*1000);
        elbow2_clr = round(elbow2_spd*1000);
        wrist2_clr = round(wrist2_spd*1000);
        hand3_clr = round(hand3_spd*1000);
        
        shoulder_clr(shoulder_clr > 200) = 200;
        elbow2_clr(elbow2_clr > 200) = 200;
        wrist2_clr(wrist2_clr > 200) = 200;
        hand3_clr(hand3_clr > 200) = 200;
        
        myColorMap = parula(250);
        figure();
        %pbaspect([2 1 1.5]);
        view(3)
        view(50,30);
        hold on
        grid on
        %plot the dots of the markers/landmarks

        %plot the link between the dots
        for i = 1:length(shoulder_pos) %iterate through each frame
            %plot each single landmark on each frame
            %scatter3(shoulder_pos(i,1),shoulder_pos(i,2),shoulder_pos(i,3),marker_size,'filled','CData',myColorMap(shoulder_clr(i),:));
            %scatter3(elbow2_pos(i,1),elbow2_pos(i,2),elbow2_pos(i,3),marker_size,'filled','CData',myColorMap(elbow2_clr(i),:));
            %scatter3(wrist2_pos(i,1),wrist2_pos(i,2),wrist2_pos(i,3),marker_size,'filled','CData',myColorMap(wrist2_clr(i),:));
            %scatter3(hand3_pos(i,1),hand3_pos(i,2),hand3_pos(i,3),marker_size,'filled','CData',myColorMap(hand3_clr(i),:));

            %plot the lines between adjacent body landmarks for each frame
            shoulder_elbow = [shoulder_pos(i,1),shoulder_pos(i,2),shoulder_pos(i,3);
                              elbow2_pos(i,1),elbow2_pos(i,2),elbow2_pos(i,3)];
            elbow_wrist = [elbow2_pos(i,1),elbow2_pos(i,2),elbow2_pos(i,3);
                           wrist2_pos(i,1),wrist2_pos(i,2),wrist2_pos(i,3)];
            wrist_hand = [wrist2_pos(i,1),wrist2_pos(i,2),wrist2_pos(i,3);
                          hand3_pos(i,1),hand3_pos(i,2),hand3_pos(i,3)];
            plot3(shoulder_elbow(:,1),shoulder_elbow(:,2),shoulder_elbow(:,3),'LineWidth',line_width,'Color','k');
            plot3(elbow_wrist(:,1),elbow_wrist(:,2),elbow_wrist(:,3),'LineWidth',line_width,'Color','k');
            plot3(wrist_hand(:,1),wrist_hand(:,2),wrist_hand(:,3),'LineWidth',line_width,'Color','k');
            
            %plot the lines (link) for each landmark between two adjacent
            %frames
            if i < length(shoulder_pos)
                adjacent_shoulder = [shoulder_pos(i,1),shoulder_pos(i,2),shoulder_pos(i,3);
                                     shoulder_pos(i+1,1),shoulder_pos(i+1,2),shoulder_pos(i+1,3)];
                adjacent_elbow = [elbow2_pos(i,1),elbow2_pos(i,2),elbow2_pos(i,3);
                                  elbow2_pos(i+1,1),elbow2_pos(i+1,2),elbow2_pos(i+1,3)];
                adjacent_wirst = [wrist2_pos(i,1),wrist2_pos(i,2),wrist2_pos(i,3);
                                  wrist2_pos(i+1,1),wrist2_pos(i+1,2),wrist2_pos(i+1,3)];
                adjacent_hand = [hand3_pos(i,1),hand3_pos(i,2),hand3_pos(i,3);
                                 hand3_pos(i+1,1),hand3_pos(i+1,2),hand3_pos(i+1,3)];
                %DASPECT
                %daspect
                %pbaspect
                %pbaspect('auto')
%                 if i_task == 1 %RT2D
%                     pbaspect([5 1 4]);
                 pbaspect([1 1 1]);
%                 else
%                     pbaspect([2 1 1.5]);
                 pbaspect([1 1 1]);
%                 end
                %pbaspect([2 1 10]);
                %arrow3(adjacent_shoulder(1,:),adjacent_shoulder(2,:),'r2',4);
                %arrow3(adjacent_elbow(1,:),adjacent_elbow(2,:),'r2',3);
                %arrow3(adjacent_wirst(1,:),adjacent_wirst(2,:),'r2',4);           
                %arrow3(adjacent_hand(1,:),adjacent_hand(2,:),'r2',1.5);   
                %plot3(adjacent_shoulder(:,1),adjacent_shoulder(:,2),adjacent_shoulder(:,3),'Color','r');
                %plot3(adjacent_elbow(:,1),adjacent_elbow(:,2),adjacent_elbow(:,3),'Color','r');
                %plot3(adjacent_wirst(:,1),adjacent_wirst(:,2),adjacent_wirst(:,3),'Color','r');           
                %plot3(adjacent_hand(:,1),adjacent_hand(:,2),adjacent_hand(:,3),'Color','r');          
            end
            
        end
        if i_task == 1
            xlim([-30,30]);
            ylim([-30,30]);
            zlim([-30,5]);
        else
            xlim([-15,20]);
            ylim([-15,20]);
            zlim([-30,5]);
        end
        
        xlabel("X Axis");
        ylabel("Y Axis");
        zlabel("Z Axis");
        
        %example of decimals in color bars
%         [c, h] = contour(peaks/30);
%         clabel(c, h)
%         cb = colorbar
%         set(cb,'YTickLabel',{'-0.20';'-0.15';'-0.10';'-0.05';' 0.00';'0.05 ';' 0.10';' 0.15';' 0.20';' 0.25';})
        %cb = colorbar;
        %cb.Ticks = 0:0.1:2.5;
        %cb.TickLabels = num2cell(string(0:0.1:2.5));
        %cb.Ticks = 0:1:8;
        %cb.TickLabels = num2cell(string(0:1:8));
        %cb.Ticks = linspace(0,3,25);
        %cb.TickLabels = num2cell(0:0.3:2.5);
        %ylabel(cb,"speed (m/s)");
        
        
        
        %hold off
        %legend show
        
        
        
        
        %start of plotting hand marker distribution, and correlation
        %between markers
        for i_marker = 1:numel(markernames) %loop for each marker
            % plot hand position for each task
            figure(f_workspace(i_marker));
            markername = markernames{i_marker};
            title(markername)
            dlc_idx = [find((strcmpi(td_list{task_idx}.dlc_pos_names,[markername,'_x']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,[markername,'_y']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,[markername,'_z'])))];
            
            temp_start = 1;
            temp_end = 4000;            
            %temp_start = 1;
            %temp_end = length(td_list{task_idx}.dlc_pos);
            
            %plot3(td_list{task_idx}.dlc_pos(:,dlc_idx(1)),td_list{task_idx}.dlc_pos(:,dlc_idx(2)),td_list{task_idx}.dlc_pos(:,dlc_idx(3)),...
            %    '.','linestyle','none','color',getColorFromList(1,i_task-1)); hold on;
            plot3(td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(1)),td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(2)),td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(3)),...
                '.','linestyle','none','color',getColorFromList(1,i_task-1)); hold on;
            %plot3(td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(1)),td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(2)),td_list{task_idx}.dlc_pos(temp_start:temp_end,dlc_idx(3)),...
            %'color',getColorFromList(1,i_task-1)); hold on;
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
        %end of plotting hand marker distribution, and correlation
        %between markers
        
        
        
        % for 2D and 3D task, compute correlation between different arm
        % markers
        figure(f_corr)
        set(gcf,'position',[0,400,2000,400])
        %corr_markernames = {'shoulder','elbow1','hand2'};
        corr_markernames = {'elbow2','hand2'};
        dlc_idx = [];
        for i_marker = 1:numel(corr_markernames)
            dlc_idx = [dlc_idx,find((strcmpi(td_list{task_idx}.dlc_pos_names,[corr_markernames{i_marker},'_x']))),...
                find((strcmpi(td_list{task_idx}.dlc_pos_names,[corr_markernames{i_marker},'_y']))),...
                ];
                %find((strcmpi(td_list{task_idx}.dlc_pos_names,[corr_markernames{i_marker},'_z'])))
        end
        dlc_data = td_list{task_idx}.dlc_pos(:,dlc_idx);
        dlc_data = dlc_data(~any(isnan(dlc_data),2),:);
        
        %Calculate velocity for each axis based on position data
        dlc_data_velocity = [diff(dlc_data)./td_list{1}.bin_size];
        %dlc_data_velocity = [diff(dlc_data./td_list{1}.bin_size)];
        %dlc_data_velocity = [diff(dlc_data)];
        %dlc_data_velocity = abs(dlc_data_velocity);
        dlc_data_speed = zeros(length(dlc_data)-1,numel(corr_markernames));
        dlc_data_diff = diff(dlc_data);
        
        %size(dlc_data_diff);
        for i_marker = 1:numel(corr_markernames)
            x = dlc_data_diff(:,i_marker*2-1);
            y = dlc_data_diff(:,i_marker*2-0);
            %x = dlc_data_diff(:,i_marker*3-2);
            %y = dlc_data_diff(:,i_marker*3-1);
            %z = dlc_data_diff(:,i_marker*3-0);
            %s = sqrt(x.^2 + y.^2 + z.^2);
            s = sqrt(x.^2 + y.^2);
            dlc_data_speed(:,i_marker) = s;
        end
        
        %corr_data = corr(dlc_data);
        corr_data = corr(dlc_data_velocity)
        %corr_data = corr(dlc_data_speed);
        corr_speed = corr(dlc_data_speed)
        
        %Save the correlation data for a scatterplot
        if task_idx == task_3d_idx
            corr_V_RT3D = corr_data;
            corr_S_RT3D = corr_speed;
        else
            corr_V_RT2D = corr_data;
            corr_S_RT2D = corr_speed;
        end
        
        subplot(1,3,i_task)
        imagesc(corr_data,[-1,1]);
        %imagesc(corr_data,[0,1]);
        %set(gca, 'xTick', [1:9]);
        
        %hardcode the axis names
        %marker_names = {'elbow X','elbow Y','elbow Z','hand X','hand Y','hand Z'};
        %set(gca,'xtick',[1:6],'xticklabel',marker_names);
        %set(gca,'ytick',[1:6],'yticklabel',marker_names);
        
        marker_names = {'elbow X','elbow Y','hand X','hand Y'};
        set(gca,'xtick',[1:4],'xticklabel',marker_names);
        set(gca,'ytick',[1:4],'yticklabel',marker_names);
        
        %marker_names = {'elbow','hand'};
        %set(gca,'xtick',[1:2],'xticklabel',marker_names);
        %set(gca,'ytick',[1:2],'yticklabel',marker_names);
        
        colormap(gray(8))
        b=colorbar();
        kin_corr(:,:,i_task) = corr_data;
        title(title_str);
    end
    
    save("kin_corr.mat","kin_corr")
    % plot correlation difference for 2D and 3D tasks, format figure
    figure(f_corr)
    subplot(1,3,3)
    corr_diff = kin_corr(:,:,1) - kin_corr(:,:,2);
    
    %saturate the correlation difference between RT3D and RT2D tasks
    corr_diff(corr_diff>0) = 1;
    corr_diff(corr_diff<0) = -1;
    
    imagesc(corr_diff,[-1,1]);
    set(gca, 'xTick', [1:4]);
    
    %das ist hardcoding
    %marker_names = {'elbow X','elbow Y','elbow Z','hand X','hand Y','hand Z'};
    %set(gca,'xtick',[1:6],'xticklabel',marker_names);
    %set(gca,'ytick',[1:6],'yticklabel',marker_names);
    
    marker_names = {'elbow X','elbow Y','hand X','hand Y'};
    set(gca,'xtick',[1:4],'xticklabel',marker_names);
    set(gca,'ytick',[1:4],'yticklabel',marker_names);
    
    %marker_names = {'elbow','hand'};
    %set(gca,'xtick',[1:2],'xticklabel',marker_names);
    %set(gca,'ytick',[1:2],'yticklabel',marker_names);
    
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
    
    
    %Make another scatter plot with the same correlation data
    %corr_data
    %corr_speed
    correlation_scatter_size = 100
    subplot();
    %scatter(corr_V_RT2D(1,2),corr_V_RT3D(1,2),correlation_scatter_size,'b','filled');
    %scatter(corr_V_RT2D(1,3),corr_V_RT3D(1,3),correlation_scatter_size,'r','filled');
    scatter(corr_V_RT2D(1,3),-0.05,correlation_scatter_size,'r','filled');
    corr_V_RT2D(1,3)
    corr_V_RT3D(1,3)
    hold on
    scatter(corr_V_RT2D(2,4),corr_V_RT3D(2,4),correlation_scatter_size,'c','filled');
    corr_V_RT2D(2,4)
    corr_V_RT3D(2,4)
    scatter(corr_S_RT2D(1,2),corr_S_RT3D(1,2),correlation_scatter_size,'k','filled');
    
    %Set the lines on the graph
    h(1) = plot([0 0],[-0.1 1],'k-','linewidth',0.5);
    h(2) = plot([-0.1 1],[0 0],'k-','linewidth',0.5);
    h(3) = plot([-0.1 1],[-0.1 1],'k--','linewidth',0.5);
    set( get( get( h(1), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    set( get( get( h(2), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    set( get( get( h(3), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    lgd = legend('Elbow-V-X - Hand-V-X',...
    'Elbow-V-Y - Hand-V-Y',...
    'Elbow-S - Hand-S',...
    '','','')
    lgd.Location = 'northwest';
    xlim([-0.1,1]);
    ylim([-0.1,1]);
    xlabel('RT2D Correlation');
    ylabel('RT3D Correlation');
    title('Correlation Comparison Between RT2D and RT3D tasks');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    %scatter(corr_speed,'r')
    %hold off
    
%         if task_idx == task_3d_idx
%             corr_V_RT3D = corr_data;
%             corr_S_RT3D = corr_speed;
%         else
%             corr_V_RT2D = corr_data;
%             corr_S_RT2D = corr_speed;
%         end
    



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