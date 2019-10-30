%% makes stim 200uA staircase data figures
    monkey_names = {'Han','Duncan'};
    folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\SingleElecs_200uA\';
    alpha_val = 0.25;
    cd(folderpath);
    use_std_err = 0; % otherwise will use standard deviations
    
    color_idx = [3,2]; %chan 19 first, then chan 42 for Duncan
    stim_amps = [40:40:200];
    offset = [-0.0005,0.0005];
    
    p_val_list = [];
%% plot stim staircase (train length, freq, amp) for each monkey (6 figures total)
    counter = 1;
    data_all = []; % min stim, bump
    num_trials_all = [];
    
    f=figure();
    f.Name = '200uA_summary_fig';
    for monk = monkey_names(1)
    % train length data
        file_list = dir([monk{1},'*.mat']);
        
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
            end
            
            num_trials_all = [num_trials_all,num_trials];
            
            min_stim_idx = numel(mean_rt);

            [h,p] = ttest2(data.cueInfo(end).rt,bump_data_rt,'vartype','unequal');
            
            p_val_list(end+1,1) = p;
            
            data_all(end+1,:) = [mean_rt(min_stim_idx),mean_bump];
            % make plots
            errorbar(mean_bump+offset(mod(file_num,2)+1),mean_rt(min_stim_idx),std_rt(min_stim_idx),std_rt(min_stim_idx),std_bump,std_bump,...
                'marker','.','color',getColorFromList(1,color_idx(counter)),'markersize',28,'linewidth',1)

        end

        counter = counter + 1;
    end

    % format figure
    unity_line = plot([0,1],[0,1],'k--','linewidth',2)
    uistack(unity_line,'bottom');
    xlim([0.1,0.28])
    ylim([0.1,0.28])
    formatForLee(gcf)
    xlabel('Bump RT (s)')
    ylabel('Stim RT (s)')
    set(gca,'fontsize',14)
    
%%

data_all
