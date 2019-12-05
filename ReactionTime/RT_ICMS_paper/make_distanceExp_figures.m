%% makes bump data figures
    monkey_names = {'Han','Duncan'};
    folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Distance_data\';
    alpha = 0.25;
    cd(folderpath);
    use_std_err = 0; % otherwise will use standard deviations
    group_size = [2,3];
    markers = {'','.','o'};
    marker_size = [0,20,12];
    
    color_idx = [3,2];

%% make scatter plot figure
    NEAR_ALL = [];
    FAR_ALL = [];

    counter = 1;
    f=figure();
    f.Name = 'DistanceExp_scatterSummary';

    for monk = monkey_names(1)
    % train length data
        file_list = dir([monk{1},'*distanceExpAll*.mat']);
        load(file_list.name);
    % collect mean data
        mean_near = [];
        mean_diff = [];
        mean_far = [];
        num_elecs = [];
        std_near = [];
        std_far = [];
        
        for group = unique(data_all.group)
            mean_near(group) = mean(data_all.rt(data_all.group==group & data_all.dist==1));
            %std_near(group) = std(data_all.rt(data_all.group==group & data_all.dist==1))/sqrt(numel(data_all.rt(data_all.group==group & data_all.dist==1)));
            std_near(group) = std(data_all.rt(data_all.group==group & data_all.dist==1));
            mean_far(group) = mean(data_all.rt(data_all.group==group & data_all.dist==0));
            %std_far(group) = std(data_all.rt(data_all.group==group & data_all.dist==0))/sqrt(numel(data_all.rt(data_all.group==group & data_all.dist==0)));
            std_far(group) = std(data_all.rt(data_all.group==group & data_all.dist==0));
            mean_diff(group) = mean_near(group)-mean_far(group);

            idx = find(data_all.group == group & data_all.dist == 1,1,'first');
            num_elecs(group) = numel(data_all.chans{idx});
            num_trials_all = [num_trials_all, sum(data_all.group==group & data_all.dist==1), sum(data_all.group==group & data_all.dist==0)];
        end
        NEAR_ALL = [NEAR_ALL, mean_near];
        FAR_ALL = [FAR_ALL, mean_far];
        % scatter plot 

        for i = group_size
            errorbar(mean_near(num_elecs == i),mean_far(num_elecs == i),...
                std_far(num_elecs == i),std_far(num_elecs == i),...
                std_near(num_elecs == i),std_near(num_elecs == i),...
                markers{i},'markersize',marker_size(i),'color',getColorFromList(1,color_idx(counter)),'linewidth',1)
            hold on
        end

    counter = counter+1;   
    end
    
    % format figure
    hold on
    temp = plot([0,0.5],[0,0.5],'k--','linewidth',1.5)
    uistack(temp,'bottom')
    xlabel('Adjacent electrodes'' RT (s)')
    ylabel('Non-adjacent electrodes'' RT(s)')
    formatForLee(gcf)
    set(gca,'fontsize',16)
    xlim([0.09,0.32])
    ylim([0.09,0.32])

    
%% make best elec vs. together plot
    TOGETHER_ALL = [];
    FASTEST_ALL = [];

    counter = 1;
    f=figure();
    f.Name = 'DistanceExp_BestVsTogether';

    for monk = monkey_names
    % train length data
        file_list_dist = dir([monk{1},'*distanceExpAll*.mat']);
        file_list_single = dir([monk{1},'*single*']);
        load(file_list_dist.name);
        load(file_list_single.name);
    % collect mean data
        mean_together = [];
        mean_fastest = [];
        num_elecs = [];
        std_together = [];
        std_fastest = [];
        
        for group = unique(data_all.group)
            mean_together(group) = mean(data_all.rt(data_all.group==group));
            std_together(group) = std(data_all.rt(data_all.group==group));

            chans = data_all.chans{find(data_all.group==group,1,'first')};
            min_rt = 10000;
            min_rt_idx = -1;
            for c = 1:numel(chans)
                if(~isempty(find(all_files_data.chan == chans(c))) && all_files_data.mean_rt(find(all_files_data.chan == chans(c))) < min_rt)
                    min_rt = all_files_data.mean_rt(find(all_files_data.chan == chans(c)));
                    min_rt_idx = c;
                end
            end
            
            all_files_data_idx = find(all_files_data.chan == chans(min_rt_idx));
            mean_fastest(group) = all_files_data.mean_rt(all_files_data_idx);
            std_fastest(group) = all_files_data.std_rt(all_files_data_idx);
            
            mean_diff(group) = mean_together(group)-mean_fastest(group);

            idx = find(data_all.group == group & data_all.dist == 1,1,'first');
            num_elecs(group) = numel(data_all.chans{idx});

        end
        FASTEST_ALL = [FASTEST_ALL,mean_fastest];
        TOGETHER_ALL = [TOGETHER_ALL,mean_together];
        
        % scatter plot 

        for i = group_size
            errorbar(mean_fastest(num_elecs == i),mean_together(num_elecs == i),...
                std_together(num_elecs == i),std_together(num_elecs == i),...
                std_fastest(num_elecs == i),std_fastest(num_elecs == i),...
                markers{i},'markersize',marker_size(i),'color',getColorFromList(1,color_idx(counter)),'linewidth',1)
            hold on
        end

    counter = counter+1;   
    end
    
    % format figure
    hold on
    temp = plot([0,0.5],[0,0.5],'k--','linewidth',1.5)
    uistack(temp,'bottom')
    xlabel('Fastest single electrode''s RT (s)')
    ylabel('Group RT(s)')
    formatForLee(gcf)
    set(gca,'fontsize',16)
    xlim([0.09,0.32])
    ylim([0.09,0.32])
    
    %% stats
    [~,p_dist] = ttest(NEAR_ALL,FAR_ALL)
    [~,p_together] = ttest(TOGETHER_ALL,FASTEST_ALL)
    
