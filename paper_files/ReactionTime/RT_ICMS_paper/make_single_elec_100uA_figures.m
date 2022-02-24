%% load in data from folder
    monkey_names = {'Han','Duncan'};
    folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\SingleElecs\';
    
    cd(folderpath);
    mean_trials = [];
    std_trials = [];
    p_vals_to_bump = [];
    p_vals_all_electrodes = [];
    lm = {};
    rotation_angles = [0,0]; %
    anova_data = [];
    anova_data.session = [];
    anova_data.ICMS = [];
    anova_data.bump = [];
    anova_data.rt = [];
    
    for monk = monkey_names(1)
        monk = monk{1};
        if(strcmpi(monk,'Duncan'))
            color_use = getColorFromList(1,2);
        else
            color_use = getColorFromList(1,3);
        end
        data_list = dir(['*',monk,'*.mat']);
        load(data_list(1).name);
    % dot plot of RTs (instead of a histogram?)
        f=figure();
        hold on;
        f.Position = [278.6000 95.4000 640 420];

        f.Name = [monk,'_singleElectrodes_dotPlot'];

        ax=gca;
        % get indexes of each day
        filename_list = unique(all_files_data.filename);
        end_idx = 0; % will be set to start idx in for loop, start at beginning
        for fnum = 1:numel(filename_list)
            % get data for this day
            start_idx = end_idx+1;
            end_idx = max(find(strcmpi(all_files_data.filename,filename_list{fnum})));
            x_data = (start_idx:1:end_idx) + fnum - 1;

            % plot bump bar for this day
            % bump data [0.1503,0.1401,0.1489,0.1450,0.1476,0.1659,0.1550,0.1639,0.1650,0.1603,0.1439] for [10/15,10/16,10/17,10/19,10/22,11/15,11/15,11/19,11/21,11/27,11/28] in Han
            % std err bump data [0.0017,0.003,0.0029,0.0032,0.0035,0.0028,0.0027,0.0022,0.0015,0.003,0.0034] for
            % [10/15,10/16,10/17,10/19,10/22,11/15,11/16,11/19,11/21,11/27] in Han

            % Duncan: 20190308: 0.1803+- 0.0102 (std dev) for 16 trials
            % 20190309: 0.1773+-0.0056,17 trials
            % 20190311: 0.178+-0.0118, 19trials
            if(strcmpi(monk,'Han'))
                bump_list = [0.1503,0.1401,0.1489,0.1450,0.1476,0.1659,0.1550,0.1639,0.1650,0.1603,0.1439]; % Han
                bump_num_trials = [10,17,17,24,12,11,27,20,16,33,18];
                bump_std_list = [0.0017,0.003,0.0029,0.0032,0.0035,0.0028,0.0027,0.0022,0.0015,0.003,0.0034].*...
                    sqrt(bump_num_trials); % Han

            elseif(strcmpi(monk,'Duncan'))
                bump_list = [0.1803,0.1773,0.178];
                bump_num_trials = [16,17,19];
                bump_std_list = [0.0102,0.0056,0.0118]; 
            end
            bump_bar = fill([x_data(1)-1,x_data(end)+1,x_data(end)+1,x_data(1)-1,x_data(1)-1],...
                bump_list(fnum)+bump_std_list(fnum).*[-1,-1,1,1,-1],...
                'k','EdgeColor','none','FaceAlpha',0.25);
            uistack(bump_bar,'bottom')
            plot([x_data(1)-1,x_data(end)+1],[bump_list(fnum),bump_list(fnum)],'k','linewidth',1.5)

            % plot visual bar for this day
            % visual data [0.2141 0.2004 0.2170 0.2281, 0.2178,0.2199,0.235,0.2122,0.247,0.2212,0.2133] for [10/15,10/16,10/17,10/19,10/22,11/15,11/16,11/19,11/21,11/27,11/28] in Han
            % std err vis data [0.0076,0.0099,0.0074,0.0089,0.0109,0.0041,0.0105,0.0043,0.0068,0.0075,0.0074] for
            % [10/15,10/16,10/17,10/19,10/22,11/15,11/16,11/19,11/21,11/27,11/28] in Han

            % duncan: 20190308: 0.4338 +- 0.0613, 18 trials
            % 20190309: 0.3925 +- 0.0439, 18 trials
            % 20190311: 0.4052 += 0.0356, 18 trials
            if(strcmpi(monk,'Han'))
                vis_list = [0.2141 0.2004 0.2170 0.2281, 0.2178,0.2199,0.235 0.2122,0.247,0.2212,0.2133]; % Han
                vis_num_trials = [14,12,18,16,17,15,15,15,20,16,10];
                vis_std_list = [0.0076,0.0099,0.0074,0.0089,0.0109,0.0041,0.0105 0.0043,0.0068,0.0075,0.0074].*...
                    sqrt(vis_num_trials); % Han
            elseif(strcmpi(monk,'Duncan'))
                vis_list = [0.4338,0.3925,0.4025];
                vis_num_trials = [18,18,18];
                vis_std_list = [0.0613,0.0439,0.0356];
            end
            vis_bar = fill([x_data(1)-1,x_data(end)+1,x_data(end)+1,x_data(1)-1,x_data(1)-1],...
            vis_list(fnum)+vis_std_list(fnum).*[-1,-1,1,1,-1],...
            'k','EdgeColor','none','FaceAlpha',0.25);
            uistack(vis_bar,'bottom')
            plot([x_data(1)-1,x_data(end)+1],[vis_list(fnum),vis_list(fnum)],'k--','linewidth',1.5)
            
            % plot data for this day
%             plot(x_data,all_files_data.mean_rt(start_idx:end_idx),'.','color',color_use,'markersize',18) 
            for i = 1:numel(x_data)
                errorbar(x_data(i),all_files_data.mean_rt(i+start_idx-1),all_files_data.std_rt(i+start_idx-1),'.','markersize',18,'color',color_use,'linewidth',1);
            end
            % plot dashed line to denote end of a day
            if(fnum < numel(filename_list)) % don't do this for the last one
                plot(x_data(end)+1+[0,0],[0,1],'--','color',[0.5,0.5,0.5],'linewidth',1.5)
            end

            % get stats for the day
            for i = 1:numel(x_data)
                rt_data = all_files_data.rt_all{i+start_idx-1};
                bump_data = all_files_data.bump_all{i+start_idx-1};
                
                [h,p] = ttest2(rt_data,bump_data,'vartype','unequal');
                p_vals_to_bump(end+1,:) = [p,mean(rt_data) < mean(bump_data)];
                
                anova_data.rt(end+1,1) = all_files_data.mean_rt(i+start_idx-1);
                anova_data.session(end+1,1) = fnum;
                anova_data.ICMS(end+1,1) = 1;
            end
            
            anova_data.rt(end+1,1) = mean(bump_data);
            anova_data.session(end+1,1) = fnum;
            anova_data.ICMS(end+1,1) = 0;
            anova_data.bump(end+1,1) = 1;
            
        end

        ylim([0,0.5])
        xlim([0,x_data(end)+1])
    %     set(gca,'XTickLabel',{});
    %     set(gca,'XTick',[]);
        xlabel('Electrode');
        ylabel('RT (s)')
        formatForLee(gcf)
        ax.FontSize = 14;
        
        
        % get mean and std of num trials
        trial_length = [all_files_data.num_trials, bump_num_trials, vis_num_trials];
%         num_trials_all = [num_trials_all, trial_length];
        mean_trials(end+1,1) = mean(trial_length);
        std_trials(end+1,1) = std(trial_length);
    
        num_trials_all = [num_trials_all,trial_length];
        
        
        
    % vertical histogram of RTs (meant to go next to above plot)
        binSize = 0.01;
        binEdges = 0.1:binSize:0.4;
        f=figure();
        f.Name = [monk,'_singleElectrodes_vertHist'];
        f.Position = [730.6000 155.4000 241.6000 420];
        ax=gca;

        % histogram
        binCounts = histcounts(all_files_data.mean_rt,binEdges)/numel(all_files_data.chan);
        h = barh(binEdges(1:end-1)+binSize/2,binCounts,'BarWidth',1);
        h.FaceAlpha = 1;
        h.FaceColor = color_use;
        h.EdgeColor = 'k';
%         ax.YAxis.TickLabels = {};

        hold on
            % plot bump bar and visual bar
        vis_mean = sum(vis_list.*vis_num_trials)/sum(vis_num_trials);
        vis_pooled_std = sqrt(sum((vis_std_list.^2).*(vis_num_trials-1))/...
            sum(vis_num_trials-1));

        vis_bar = fill([0,1,1,0,0],...
            vis_mean+vis_pooled_std.*[-1,-1,1,1,-1],...
            'k','EdgeColor','none','FaceAlpha',0.3);
        uistack(vis_bar,'bottom')
        plot([0,1],[vis_mean,vis_mean],'k--','linewidth',1.5)

        bump_mean = sum(bump_list.*bump_num_trials)/sum(bump_num_trials);
        bump_pooled_std = sqrt(sum((bump_std_list.^2).*(bump_num_trials-1))/...
            sum(bump_num_trials-1));

        % do all electrodes stats
        % paired test, find difference between electrode and corresponding
        % bump
        stat_rt = [];
        stat_bump = [];
        for i = 1:numel(all_files_data.mean_rt)
            stat_rt(end+1) = all_files_data.mean_rt(i); 
            stat_bump(end+1) = mean(all_files_data.bump_all{i});
        end
        stat_bump = unique(stat_bump); % only use 1 bump per session and do not pair data
        
        [h,p] = ttest2(stat_rt,stat_bump,'vartype','unequal');
        p_vals_all_electrodes(end+1) = p;
        
        
        bump_bar = fill([0,1,1,0,0],...
                bump_mean+bump_pooled_std.*[-1,-1,1,1,-1],...
                'k','EdgeColor','none','FaceAlpha',0.3);
        uistack(bump_bar,'bottom')
        plot([0,1],[bump_mean,bump_mean],'k','linewidth',1.5) 


        xlabel('Proportion of electrodes');
        ylabel('RT (s)');
        formatForLee(gcf)
        set(gca,'fontsize',14)
        ylim([0.,0.5])
        xlim([0,0.2])
        
% plot heatmap of RT across array
        pos_all = [];
        if(strcmpi(monk,'Han'))
            mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
            bad_elecs = [];
        elseif(strcmpi(monk,'Duncan'))
            mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
            bad_elecs = [63,96];
        end
        
        map_data = loadMapFile(mapFileName);
        c_map = inferno();
        c_map = flip(c_map,1);
        rt_min = 0.1;
        rt_max = 0.4;
        
        
        f=figure();
        f.Name = [monk,'_singleElectrodes_arrayMap'];
        hold on
        for chan_idx = 1:numel(all_files_data.chan)
            map_idx = find(all_files_data.chan(chan_idx) == map_data.chan);
            mean_rt = all_files_data.mean_rt(chan_idx);
            c_map_idx = ceil((mean_rt-rt_min)/(rt_max-rt_min)*size(c_map,1));
            rt_color = c_map(min(max(c_map_idx,1),size(c_map,1)),:);
            pos_all(chan_idx,:) = [map_data.row(map_idx),11-map_data.col(map_idx)];
            rectangle('Position',[map_data.row(map_idx),11-map_data.col(map_idx),1,1],'FaceColor',rt_color, 'EdgeColor','none');
            hold on
        end
        plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5)

        % plot channels with no response
        for i = 1:numel(all_files_data.no_response)
            map_idx = find(all_files_data.no_response(i) == map_data.chan);
            pos = [map_data.row(map_idx),11-map_data.col(map_idx)];
            plot([pos(1),pos(1)+1],[pos(2),pos(2)+1],'k-','linewidth',2)
            plot([pos(1),pos(1)+1],[pos(2)+1,pos(2)],'k-','linewidth',2)
        end

        % plot channels that made animal vocalize
        for i = 1:numel(bad_elecs)
            map_idx = find(bad_elecs(i) == map_data.chan);
            pos = [map_data.row(map_idx),11-map_data.col(map_idx)];
            plot([pos(1),pos(1)+1],[pos(2),pos(2)+1],'r-','linewidth',2)
            plot([pos(1),pos(1)+1],[pos(2)+1,pos(2)],'r-','linewidth',2)
        end
        
        f=gcf;
        set(gca,'visible','off')
        xlim([1,11])
        ylim([1,11])
        axis square
        b=colorbar();
        colormap(c_map);
        b.FontSize = 14;
        b.Ticks = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]; 
        b.TickDirection = 'out';
        b.TickLabels = cell(1,numel(b.Ticks));
        for i = 1:2:numel(b.Ticks)
            if(i==numel(b.Ticks))
                b.TickLabels{i,1} = strcat(num2str((i-1)*(rt_max-rt_min)/(numel(b.Ticks)-1) + rt_min));
            elseif(i==1)
                b.TickLabels{i,1} = strcat(num2str((i-1)*(rt_max-rt_min)/(numel(b.Ticks)-1) + rt_min));
            else
                b.TickLabels{i,1} = num2str((i-1)*(rt_max-rt_min)/(numel(b.Ticks)-1) + rt_min);
            end

        end
        b.Label.String = 'RT (s)';
        b.Label.FontSize = 16;
        
        % linear model for each monkey
        if(strcmpi(monk,'Han')==1)
            rot = rotation_angles(1);
        else
            rot = rotation_angles(2);
        end
        pos_all = [pos_all(:,1)*cos(rot) + pos_all(:,2)*sin(rot), -pos_all(:,1)*sin(rot)+pos_all(:,2)*cos(rot)];
        
        lm{end+1} = fitlm(pos_all,all_files_data.mean_rt,'y~x1+x2');
        
        
    end
