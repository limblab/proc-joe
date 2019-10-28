%% makes 16 elec amp staircase data figures
    monkey_names = {'Han','Duncan'};
    folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Multielec_staircase\';
    cd(folderpath);
    use_std_err = 0; % otherwise will use standard deviations
    
    color_idx = [3,2];
    markers = {'.','x','s','o'};
    marker_size = [28,12,12,12];
    USED_STIM_AMP = [30,35,40,50];

    offset = [-0.0005,0.0005];
    p_vals = [];
    mean_difference = [];
%% plot stim staircase (train length, freq, amp) for each monkey (6 figures total)
    counter = 1;
    data_all = []; % min stim, bump
    f=figure();
    f.Name = '16elecs_summary_figure';
   
%     num_trials_all = [];
    
    for monk = monkey_names(2)
        stat_bump = [];
        stat_stim = [];
        
        file_list = dir([monk{1},'*data.mat']);
        
        hold on
        for file_num = 1:numel(file_list)
            load(file_list(file_num).name);
            
            % bump cue data
            bump_data_rt = [];
            if(any([data.cueInfo.bumpMag] > 0 & [data.cueInfo.stimCode] == -1))
                bump_data_rt = data.cueInfo([data.cueInfo.bumpMag] > 0 & [data.cueInfo.stimCode] == -1).rt;
            end
            
            data.cueInfo([data.cueInfo.bumpMag] > 0) = [];
            data.cueInfo([data.cueInfo.stimCode] < 0) = [];
            
            mean_bump = mean(bump_data_rt);
            std_bump = std(bump_data_rt);
            bump_num_trials = numel(bump_data_rt);
            if(use_std_err)
                std_bump = std_bump/sqrt(bump_num_trials);
            end
            
            num_trials_all = [num_trials_all,bump_num_trials];
            
            % data holds the bump data with bump mag as a variable
            mean_rt = []; std_rt = []; num_trials = []; stim_amps_all = []; rt_data_all = [];
            for d = 1:numel(data.cueInfo)
                mean_rt(d) = mean(data.cueInfo(d).rt);
                std_rt(d) = std(data.cueInfo(d).rt);
                num_trials(d) = numel(data.cueInfo(d).rt);
                if(use_std_err)
                    std_rt(d) = std_rt(d)/sqrt(num_trials(d));
                end
                stim_amps_all = [stim_amps_all,STIM_PARAMS(d)*ones(1,numel(data.cueInfo(d).rt))];
                rt_data_all = [rt_data_all,data.cueInfo(d).rt];
            end

            num_trials_all = [num_trials_all,num_trials];
            
            [~,min_stim_idx] = min(mean_rt);
            stim_amp = STIM_PARAMS(min_stim_idx);
            marker_idx = find(stim_amp == USED_STIM_AMP);
            
            stat_bump(end+1) = mean_bump;
            stat_stim(end+1) = mean_rt(min_stim_idx);
            
            data_all(end+1,:) = [mean_rt(min_stim_idx),mean_bump];
            % make plots
            errorbar(mean_bump+offset(mod(file_num,2)+1),mean_rt(min_stim_idx),std_rt(min_stim_idx),std_rt(min_stim_idx),std_bump,std_bump,...
                'marker',markers{marker_idx},'color',getColorFromList(1,color_idx(counter)),'markersize',marker_size(marker_idx),'linewidth',1)

        end

        [~,p_vals(end+1)] = ttest(stat_stim,stat_bump);
        mean_difference(end+1) = mean(stat_stim-stat_bump);
        counter = counter + 1;
    end

    % format figure
    unity_line = plot([0,1],[0,1],'k--','linewidth',2);
    uistack(unity_line,'bottom');    
    xlim([0.1,0.25])
    ylim([0.1,0.25])
    formatForLee(gcf)
    xlabel('Bump RT (s)')
    ylabel('Stim RT (s)')
    set(gca,'fontsize',14)
    
