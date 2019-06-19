%% makes bump data figures
    monkey_names = {'Han','Duncan'};
    folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\SingleElecs_staircase\';
    alpha_val = 0.25;
    cd(folderpath);
    use_std_err = 0; % otherwise will use standard deviations
    
    color_idx = [0,1;0,1]; % chan 35 = 0, chan 37 = 1 for Duncan, chan 20 = 0 , chan 11 = 1,  chan 62 = 4 for Han
    amp_offset = [-0.75,0.75];
%% plot stim staircase for train length data (train length, freq, amp) for each monkey (6 figures total)
    amps = {[10:10:100]};
    
    ANOVA_data = [];
    ANOVA_data.rt = [];
    ANOVA_data.amp = [];
    ANOVA_data.chan_monk = [];
    
    counter = 1;
    r2_all = [];
    
    for monk = monkey_names
        % train length data
        file_list = dir([monk{1},'*EXAMPLE*amp*']);
        
        f=figure();
        f.Name = [monk{1},'_singleElec_ampStaircase'];
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
            
            % remove stim codes that weren't felt
            amp_idx = find(cellfun(@numel,amps) == numel(data.cueInfo));
            amps_use = amps{amp_idx};
            
            data.cueInfo(isnan([data.cueInfo.percent_respond])) = [];
            % data holds the bump data with bump mag as a variable
            
            mean_rt = []; std_rt = []; num_trials = []; amp_all = []; rt_data_all = []; a = [];
            for d = 1:numel(data.cueInfo)
                mean_rt(d) = mean(data.cueInfo(d).rt);
                std_rt(d) = std(data.cueInfo(d).rt);
                num_trials(d) = numel(data.cueInfo(d).rt);
                a(d) = amps_use(data.cueInfo(d).stimCode+1);
                if(use_std_err)
                    std_rt(d) = std_rt(d)/num_trials(d);
                end
                amp_all = [amp_all,a(d)*ones(1,numel(data.cueInfo(d).rt))];
                rt_data_all = [rt_data_all,data.cueInfo(d).rt];
            end
            
            chan_all = (counter-1)*10 + file_num + zeros(size(a));
            monk_all = counter + zeros(size(a));
            ANOVA_data.rt = [ANOVA_data.rt,mean_rt];
            ANOVA_data.amp = [ANOVA_data.amp, a];
            ANOVA_data.chan_monk = [ANOVA_data.chan_monk, chan_all];
            
            
            % fit data
            [fitObj,gof] = fit(amp_all',rt_data_all','a*exp(b*x)+c','startPoint',[0,0,0]);
            r2_all(counter,file_num) = gof.rsquare;
            
            x_data_fit = linspace(min(amps_use*0.9),max(amps_use)*1.1,100);
            y_data_fit = fitObj.a*exp(fitObj.b*x_data_fit)+fitObj.c;
            plot(x_data_fit,y_data_fit,'--','linewidth',2,'color',getColorFromList(1,color_idx(counter,file_num)));
            
            % make plots
            % visual bar and horizontal line
            visual_bar = fill([0,max(amps_use)*2,max(amps_use)*2,0],...
                mean_vis(1)+std_vis*[-1,-1,1,1],'k','linestyle','none');
            alpha(alpha_val);
            visual_line=plot([0,max(amps_use)*2],mean_vis+[0,0],'--','linewidth',2,'color','k');
            uistack(visual_line,'bottom');
            uistack(visual_bar,'bottom');
                        % bump bar and horz line
            if(~isempty(bump_data_rt))
                bump_bar = fill([0,max(amps_use)*2,max(amps_use)*2,0],...
                    mean_bump(1)+std_bump*[-1,-1,1,1],'k','linestyle','none');
                alpha(alpha_val);
                bump_line=plot([0,max(amps_use)*2],mean_bump+[0,0],'-','linewidth',2,'color','k');
                uistack(bump_line,'bottom');
                uistack(bump_bar,'bottom');
            end
            

            
            
            errorbar(a+amp_offset(file_num),mean_rt,std_rt,'.','color',getColorFromList(1,color_idx(counter,file_num)),...
                'markersize',28,'linewidth',2)
            hold on
            for d = 1:numel(data.cueInfo)
                dots_stim = scatter(a(d)+amp_offset(file_num)+zeros(1,numel(data.cueInfo(d).rt)),data.cueInfo(d).rt,...
                    28,'markeredgecolor','none','markerfacecolor',getColorFromList(1,color_idx(counter,file_num)))
                alpha(dots_stim,0.5)
            end
        end
        % format figure
        xlim([0,max(amps_use)*1.1]);
        ylim([0,0.5]);
        formatForLee(gcf);
        xlabel('Amplitude (\muA)');
        ylabel('RT (s)');
        
        set(gca,'fontsize',14)
        ax = gca;
        ax.XTick = [0:20:100];
        ax.XAxis.MinorTickValues = [0:10:100];

        counter = counter + 1;
    end
    
    ANOVA_data.chan_monk = categorical(ANOVA_data.chan_monk);
    
    ANOVA_data.chan_monk = ANOVA_data.chan_monk';
    ANOVA_data.rt = ANOVA_data.rt';
    ANOVA_data.amp = ANOVA_data.amp';
    %%
    keep_mask = ~isnan(ANOVA_data.rt);
    for field_name = fieldnames(ANOVA_data)'
        ANOVA_data.(field_name{1}) = ANOVA_data.(field_name{1})(keep_mask);
        
    end
   %% do ANOVA
    lm_table = table(ANOVA_data.rt,ANOVA_data.chan_monk,ANOVA_data.amp,...
        'VariableNames',{'RT','chan_monk','amp'});
    
    lm = fitlm(lm_table,'RT~chan_monk+amp')
    
    