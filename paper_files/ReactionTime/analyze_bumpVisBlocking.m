    td_vis = td_reward([td_reward.isVisualTrial]==1);
    td_bump = td_reward([td_reward.isBumpTrial]==1);


    rt_vis = ([td_vis.idx_movement_on] - [td_vis.idx_goCueTime])*td_reward(1).bin_size;
    rt_bump = ([td_bump.idx_movement_on] - [td_bump.idx_goCueTime])*td_reward(1).bin_size;

%% 135:end, 1:79 for joe; 1:54, 95:end for jul; % Han 10/09: 1:58 alone
% Han 10/08 1:67 alone, 68:119 tog, 120:end alone
% Cracle 11/09 : 1:60alone, 61:111 together, 112:end alone
% Crackle 11/10 : 11:67alone, 68:136together, 137:end
    rt_vis_alone_1 = rt_vis(11:67); 
    rt_vis_tog = rt_vis(68:136);
    rt_vis_alone_2 = rt_vis(137:end);
    
    rt_vis_block(11:67) = 1;
    rt_vis_block(68:136) = 2;
    rt_vis_block(137:numel(rt_vis)) = 3;

    
%     100:end for jul; Han 10/09: 1:64 alone
%     Han 10/08: 1:47 tog, 48:end alone Bump
%       Crackle 11/09 : 1:47 together, 48:end alone
%   Crackle 11/10 : 13:53alone, 54:92together, 93:end
    rt_bump_alone_1 = rt_bump(13:53); % 128:end, 1:66 for joe; 1:55,
    rt_bump_tog = rt_bump(54:92);
    rt_bump_alone_2 = rt_bump(93:end);

    rt_bump_block(13:53) = 1;
    rt_bump_block(54:92) = 2;
    rt_bump_block(93:numel(rt_bump)) = 3;

%% plot stuff

    vis_color = getColorFromList(1,0);
    bump_color = getColorFromList(1,1);
    marker_list = {'.','s','o'};
    marker_sizes = [8,6,4];

    figure(); hold on
    for i_mode = 1:3
        plot([td_vis(rt_vis_block==i_mode).trial_id],rt_vis(rt_vis_block==i_mode),...
            marker_list{i_mode},'color',vis_color,'markersize',marker_sizes(i_mode))
        plot([td_bump(rt_bump_block==i_mode).trial_id],rt_bump(rt_bump_block==i_mode),...
            marker_list{i_mode},'color',bump_color,'markersize',marker_sizes(i_mode))
    end

    l=legend('Visual','Bump');
    xlabel('Trial');
    ylabel('RT (s)');
    formatForLee(gcf)
    set(gca,'fontsize',14)


    figure(); hold on
    
    errorbar(-0.05,mean(rt_vis_alone_1,'omitnan'),std(rt_vis_alone_1,'omitnan'),...
        'color',vis_color,'linewidth',1.5,'marker',marker_list{1},'markersize',40)
    errorbar(0.0,mean(rt_vis_tog,'omitnan'),std(rt_vis_tog,'omitnan'),...
        'color',vis_color,'linewidth',0.5,'marker',marker_list{2},'markersize',12)
    
    if(exist('rt_vis_alone_2')>0)
        errorbar(0.05,mean(rt_vis_alone_2,'omitnan'),std(rt_vis_alone_2,'omitnan'),...
            'color',vis_color,'linewidth',1.5,'marker',marker_list{3},'markersize',12)
    end
    
    if(exist('rt_bump_alone_1')>0)
        errorbar(0.95,mean(rt_bump_alone_1,'omitnan'),std(rt_bump_alone_1,'omitnan'),...
            'color',bump_color,'linewidth',1.5,'marker',marker_list{1},'markersize',40)
    end
    errorbar(1.0,mean(rt_bump_tog,'omitnan'),std(rt_bump_tog,'omitnan'),...
        'color',bump_color,'linewidth',0.5,'marker',marker_list{2},'markersize',12)
    
    if(exist('rt_bump_alone_2')>0)
        errorbar(1.05,mean(rt_bump_alone_2,'omitnan'),std(rt_bump_alone_2,'omitnan'),...
            'color',bump_color,'linewidth',1.5,'marker',marker_list{3},'markersize',12)
    end
    
    ax=gca;
    xlabel('Cue Type');
    ylabel('RT (s)');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    ax.XTick = [0,1];
    ax.XMinorTick = 'off';
    ax.XTickLabel = {'Visual','Bump'};