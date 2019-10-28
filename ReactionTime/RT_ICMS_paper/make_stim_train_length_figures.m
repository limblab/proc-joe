%% makes bump data figures
    monkey_names = {'Han','Duncan'};
    folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\SingleElecs_staircase\';
    alpha_val = 0.25;
    cd(folderpath);
    use_std_err = 0; % otherwise will use standard deviations
    
    color_idx = [1,4;0,1]; % chan 35 = 0, chan 37 = 1 for Duncan, chan 20 = 0 , chan 11 = 1,  chan 62 = 4 for Han
    train_length_offset = [-1.5,1.5];
%% plot stim staircase for train length data (train length, freq, amp) for each monkey (6 figures total)
    train_lengths = {[25,50,75,100,125,150,200,250,300],...
        [25:25:200,250,300]};
    
    ANOVA_data = [];
    ANOVA_data.rt = [];
    ANOVA_data.tl = [];
    ANOVA_data.monk = [];
    ANOVA_data.chan = [];
    
    counter = 1;
    r2_all = [];
    
    for monk = monkey_names(2)
        % train length data
        file_list = dir([monk{1},'*EXAMPLE*trainLength*']);
        
        f=figure();
        f.Name = [monk{1},'_singleElec_trainLengthStaircase'];
        hold on
        for file_num = 1:numel(file_list)
            load(file_list(file_num).name);
            
            % bump cue data
            bump_data_rt = [];
            if(any([data.cueInfo.bumpMag] > 0 & [data.cueInfo.stimCode] == -1))
                bump_data_rt = data.cueInfo([data.cueInfo.bumpMag] > 0 & [data.cueInfo.stimCode] == -1).rt;
            end
            data.cueInfo([data.cueInfo.bumpMag] > 0) = [];
           
            mean_bump = mean(bump_data_rt);
            std_bump = std(bump_data_rt);
            bump_num_trials = numel(bump_data_rt);
            if(use_std_err)
                std_bump = std_bump/sqrt(bump_num_trials);
            end
            num_trials_all = [num_trials_all,bump_num_trials];

            % visual cue data
            vis_data_rt = [];
            if(any(isEqual([data.cueInfo.bumpMag],0) & [data.cueInfo.stimCode] == -1))
                vis_data_rt = data.cueInfo(isEqual([data.cueInfo.bumpMag],0) & [data.cueInfo.stimCode] == -1).rt;
            end
            data.cueInfo([data.cueInfo.stimCode]==-1) = [];
            
            mean_vis = mean(vis_data_rt);
            std_vis = std(vis_data_rt);
            vis_num_trials = numel(vis_data_rt);
            if(use_std_err)
                std_vis = std_vis/sqrt(vis_num_trials);
            end
            
            num_trials_all = [num_trials_all,vis_num_trials];

            
            % remove stim codes that weren't felt
            tl_idx = find(cellfun(@numel,train_lengths) == numel(data.cueInfo));
            train_lengths_use = train_lengths{tl_idx};
            
            data.cueInfo(isnan([data.cueInfo.percent_respond])) = [];
            % data holds the bump data with bump mag as a variable
            
            mean_rt = []; std_rt = []; num_trials = []; tl_all = []; rt_data_all = []; tl = [];
            for d = 1:numel(data.cueInfo)
                mean_rt(d) = mean(data.cueInfo(d).rt);
                std_rt(d) = std(data.cueInfo(d).rt);
                num_trials(d) = numel(data.cueInfo(d).rt);
                tl(d) = train_lengths_use(data.cueInfo(d).stimCode+1);
                if(use_std_err)
                    std_rt(d) = std_rt(d)/num_trials(d);
                end
                tl_all = [tl_all,tl(d)*ones(1,numel(data.cueInfo(d).rt))];
                rt_data_all = [rt_data_all,data.cueInfo(d).rt];
            end
            
            num_trials_all = [num_trials_all,num_trials];

            
            [fitObj,gof] = fit(tl_all',rt_data_all','a*exp(b*x)+c','startPoint',[0,0,0]);
            r2_all(counter,file_num) = gof.rsquare;
            
            x_data_fit = linspace(min(train_lengths_use*0.9),max(train_lengths_use)*1.1,100);
            y_data_fit = fitObj.a*exp(fitObj.b*x_data_fit)+fitObj.c;
            plot(x_data_fit,y_data_fit,'--','linewidth',2,'color',getColorFromList(1,color_idx(counter,file_num)));
            
            chan_all = (counter-1)+file_num + zeros(size(tl));
            monk_all = counter + zeros(size(tl));
            ANOVA_data.rt = [ANOVA_data.rt,mean_rt];
            ANOVA_data.tl = [ANOVA_data.tl, tl];
            ANOVA_data.monk = [ANOVA_data.monk, monk_all];
            ANOVA_data.chan = [ANOVA_data.chan, chan_all];
            
            % make plots
            % visual bar and horizontal line
            if(file_num==2)
                visual_bar = fill([0,max(train_lengths_use)*2,max(train_lengths_use)*2,0],...
                    mean_vis(1)+std_vis*[-1,-1,1,1],'k','linestyle','none');
                alpha(alpha_val);
                visual_line=plot([0,max(train_lengths_use)*2],mean_vis+[0,0],'--','linewidth',2,'color','k');
                uistack(visual_line,'bottom');
                uistack(visual_bar,'bottom');
                            % bump bar and horz line
                if(~isempty(bump_data_rt))
                    bump_bar = fill([0,max(train_lengths_use)*2,max(train_lengths_use)*2,0],...
                        mean_bump(1)+std_bump*[-1,-1,1,1],'k','linestyle','none');
                    alpha(alpha_val);
                    bump_line=plot([0,max(train_lengths_use)*2],mean_bump+[0,0],'-','linewidth',2,'color','k');
                    uistack(bump_line,'bottom');
                    uistack(bump_bar,'bottom');
                end
            end
            

            
            
            errorbar(tl+train_length_offset(file_num),mean_rt,std_rt,'.','color',getColorFromList(1,color_idx(counter,file_num)),...
                'markersize',28,'linewidth',2)
            hold on
            for d = 1:numel(data.cueInfo)
                dots_stim = scatter(tl(d)+train_length_offset(file_num)+zeros(1,numel(data.cueInfo(d).rt)),data.cueInfo(d).rt,...
                    28,'markeredgecolor','none','markerfacecolor',getColorFromList(1,color_idx(counter,file_num)));
                alpha(dots_stim,0.5)
            end
        end
        % format figure
        xlim([0,max(train_lengths_use)*1.1]);
        ylim([0,0.5]);
        formatForLee(gcf);
        xlabel('Train length (ms)');
        ylabel('RT (s)');
        
        set(gca,'fontsize',14)
        ax = gca;
        ax.XTick = [0:50:300];
        ax.XAxis.MinorTickValues = [0:25:300];

        counter = counter + 1;
    end
    
    
    ANOVA_data.chan = categorical(ANOVA_data.chan);
    ANOVA_data.monk = categorical(ANOVA_data.monk);
    
    ANOVA_data.chan = ANOVA_data.chan';
    ANOVA_data.rt = ANOVA_data.rt';
    ANOVA_data.monk = ANOVA_data.monk';
    ANOVA_data.tl = ANOVA_data.tl';
   %% do ANOVA
    lm_table = table(ANOVA_data.rt,ANOVA_data.chan,ANOVA_data.monk,ANOVA_data.tl,...
        'VariableNames',{'RT','chan','monk','tl'});
    
    lm = fitlm(lm_table,'RT~chan+monk+tl')
    