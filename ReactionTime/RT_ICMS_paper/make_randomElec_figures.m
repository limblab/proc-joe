%% makes 16 elec amp staircase data figures
    monkey_names = {'Han','Duncan'};
    folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\RandomElecs\';
    cd(folderpath);
    use_std_err = 0; % otherwise will use standard deviations
    alpha_val = 0.25;
    color_idx = repmat([1,4,0],1,5);
    offset = repmat([-0.35,0,0.35],1,5);
    num_elecs = [4,4,4,6,6,6,8,8,8,12,12,12,24,24,24];
    total_amp = [240,360,480,240,360,480,240,360,480,240,360,480,240,360,480];
%% plot stim staircase (train length, freq, amp) for each monkey (6 figures total)
    DATA_ALL = [];

    counter = 1;
    data_all = []; % min stim, bump
    for monk = monkey_names
    % train length data
        file_list = dir([monk{1},'*.mat']);
        
        f=figure();
        f.Name = [monk{1},'_randomElecs'];
        hold on
        load(file_list.name);
            
        % bump cue data
        bump_data_rt = [];
        if(any([data.cueInfo.bumpMag] > 0 & [data.cueInfo.stimCode] == -1))
            bump_data_rt = data.cueInfo([data.cueInfo.bumpMag] > 0 & [data.cueInfo.stimCode] == -1).rt;
        end

        data.cueInfo([data.cueInfo.bumpMag] > 0) = [];

        mean_bump = mean(BUMP_DATA);
        std_bump = POOLED_STANDARD_DEV_BUMP;
        if(use_std_err)
            std_bump = std_bump*sqrt(sum(1/bump_num_trials));
        end

        % visual cue data
        vis_data_rt = [];
        if(any(isEqual([data.cueInfo.bumpMag],0) & [data.cueInfo.stimCode] == -1))
            vis_data_rt = data.cueInfo(isEqual([data.cueInfo.bumpMag],0) & [data.cueInfo.stimCode] == -1).rt;
        end
        data.cueInfo([data.cueInfo.stimCode]==-1) = [];

        mean_vis = mean(vis_data_rt);
        std_vis = POOLED_STANDARD_DEV_VIS;
        if(use_std_err)
            std_vis = std_vis/sqrt(vis_num_trials);
        end
        
        % data holds the bump data with bump mag as a variable
        mean_rt = []; std_rt = []; num_trials = []; num_elecs_all = []; rt_data_all = [];
        for d = 1:numel(data.cueInfo)
            mean_rt(d) = mean(data.cueInfo(d).rt);
            std_rt(d) = std(data.cueInfo(d).rt);
            num_trials(d) = numel(data.cueInfo(d).rt);
            if(use_std_err)
                std_rt(d) = std_rt(d)/sqrt(num_trials(d));
            end
            num_elecs_all = [num_elecs_all,total_amp(d)*ones(1,numel(data.cueInfo(d).rt))];
            rt_data_all = [rt_data_all,data.cueInfo(d).rt];
%             DATA_ALL = [DATA_ALL; data.cueInfo(d).rt',num_elecs(d)+zeros(numel(data.cueInfo(d).rt),1),...
%                 total_amp(d)+zeros(numel(data.cueInfo(d).rt),1), counter+zeros(numel(data.cueInfo(d).rt),1)];
            DATA_ALL = [DATA_ALL; mean_rt(d), num_elecs(d), total_amp(d), counter];
        end
        [~,min_stim_idx] = min(mean_rt);

        % fit line to data to get significance?
%         [fitObj,gof] = fit(num_elecs_all',rt_data_all','a*exp(b*x)+c','startPoint',[0,0,0]);
%         r2_all(counter) = gof.rsquare;
% 
%         x_data_fit = linspace(min(num_elecs*0.9),max(num_elecs)*1.1,100);
%         y_data_fit = fitObj.a*exp(fitObj.b*x_data_fit)+fitObj.c;
%         plot(x_data_fit,y_data_fit,'--','linewidth',2,'color',getColorFromList(1,color_idx(f)));


        % make plots
        % visual bar and horizontal line
        if(~isempty(vis_data_rt))
            visual_bar = fill([0,max(num_elecs)*2,max(num_elecs)*2,0],...
                mean_vis(1)+std_vis*[-1,-1,1,1],'k','linestyle','none');
            alpha(alpha_val);
            visual_line=plot([0,max(num_elecs)*2],mean_vis+[0,0],'--','linewidth',2,'color','k');
            uistack(visual_line,'bottom');
            uistack(visual_bar,'bottom');
        end
            
        % bump bar and horz line
        if(~isempty(bump_data_rt))
            bump_bar = fill([0,max(num_elecs)*2,max(num_elecs)*2,0],...
                mean_bump(1)+std_bump*[-1,-1,1,1],'k','linestyle','none');
            alpha(alpha_val);
            bump_line=plot([0,max(num_elecs)*2],mean_bump+[0,0],'-','linewidth',2,'color','k');
            uistack(bump_line,'bottom');
            uistack(bump_bar,'bottom');
        end

        for d = 1:numel(data.cueInfo)
            errorbar(num_elecs(d)+offset(d),mean_rt(d),std_rt(d),...
                '.','color',getColorFromList(1,color_idx(d)),'markersize',28,'linewidth',1.5)
            hold on
            dots_stim = scatter(num_elecs(d)+offset(d)+zeros(1,numel(data.cueInfo(d).rt)),data.cueInfo(d).rt,...
                28,'markeredgecolor','none','markerfacecolor',getColorFromList(1,color_idx(d)));
            alpha(dots_stim,0.5)
        end
        
        counter = counter + 1;
        
        % format figure
        xlim([0,27])
        ylim([0,0.5])
        xlabel('Num elecs')
        ylabel('RT (s)')
        formatForLee(gcf)
        set(gca,'fontsize',14)
        ax = gca;
        
        ax.XTick = unique(num_elecs);
        ax.XAxis.MinorTickValues = [0:1:28];
    end

    %% lm with DATA_ALL
    data_table = array2table(DATA_ALL,'VariableNames',{'rt','num_elecs','total_amp','monk'});
    data_table.monk = categorical(data_table.monk);
    lm = fitlm(data_table,'rt~num_elecs+total_amp+monk')