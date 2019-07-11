    data_all.rt = [];
    data_all.group = [];
    data_all.dist = [];
    data_all.chans = {};
    group_idx = 1;

%%
    for c = 1:10
        data_all.rt = [data_all.rt,data.cueInfo(c).rt];
        data_all.group = [data_all.group,group_idx+zeros(size(data.cueInfo(c).rt))];
        data_all.dist = [data_all.dist,mod(c,2)+zeros(size(data.cueInfo(c).rt))];
        data_all.chans = [data_all.chans,repmat(electrodeList(c),1,numel(data.cueInfo(c).rt))];
        if(mod(c,2) == 0)
            group_idx = group_idx+1;
        end
    end

%%
    mean_near = [];
    mean_diff = [];
    mean_far = [];
    num_elecs = [];
    for group = unique(data_all.group)
        mean_near(group) = mean(data_all.rt(data_all.group==group & data_all.dist==1));
        std_err_near(group) = std(data_all.rt(data_all.group==group & data_all.dist==1))/sqrt(numel(data_all.rt(data_all.group==group & data_all.dist==1)));
        mean_far(group) = mean(data_all.rt(data_all.group==group & data_all.dist==0));
        std_err_far(group) = std(data_all.rt(data_all.group==group & data_all.dist==0))/sqrt(numel(data_all.rt(data_all.group==group & data_all.dist==0)));
        mean_diff(group) = mean_near(group)-mean_far(group);
        
        idx = find(data_all.group == group & data_all.dist == 1,1,'first');
        num_elecs(group) = numel(data_all.chans{idx});

    end

%% scatter plot 
    f=figure();
    f.Name = 'Duncan_distanceExp_scatterSummary';
    markers = {'','.','o'};
    marker_size = [0,20,10];
    for i = 2:3
        errorbar(mean_near(num_elecs == i),mean_far(num_elecs == i),std_err_near(num_elecs == i),...
            'horizontal',markers{i},'markersize',marker_size(i),'color','k')
        hold on
        errorbar(mean_near(num_elecs == i),mean_far(num_elecs == i),std_err_far(num_elecs == i),...
            'vertical',markers{i},'markersize',marker_size(i),'color','k')
    end
    
    hold on
    temp = plot([0,0.3],[0,0.3],'r--','linewidth',1.5)
    uistack(temp,'bottom')
    xlabel('Adjacent electrodes'' RT (s)')
    ylabel('Non-adjacent electrodes'' RT(s)')
    formatForLee(gcf)
    set(gca,'fontsize',16)
    xlim([0.1,0.25])
    ylim([0.1,0.25])
