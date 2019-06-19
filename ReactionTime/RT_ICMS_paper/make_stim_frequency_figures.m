%% makes bump data figures
    monkey_names = {'Han','Duncan'};
    folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\SingleElecs_staircase\';
    alpha_val = 0.25;
    cd(folderpath);
    use_std_err = 0; % otherwise will use standard deviations
    
    color_idx = [0,1;0,1]; % chan 35 = 0, chan 37 = 1 for Duncan, chan 20 = 0 , chan 11 = 1,  chan 62 = 4 for Han
    frequency_offset = [-2.5,2.5];
%% plot stim staircase for train length data (train length, freq, amp) for each monkey (6 figures total)
    frequencies = {[50:50:500],...
        [50:50:400]};
    
    
    ANOVA_data = [];
    ANOVA_data.rt = [];
    ANOVA_data.freq = [];
    ANOVA_data.monk = [];
    ANOVA_data.chan = [];
    
    
    counter = 1;
    r2_all = [];
    
    for monk = monkey_names
        % train length data
        file_list = dir([monk{1},'*EXAMPLE*frequency*']);
        
        f=figure();
        f.Name = [monk{1},'_singleElec_frequencyStaircase'];
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
            freq_idx = find(cellfun(@numel,frequencies) == numel(data.cueInfo));
            frequencies_use = frequencies{freq_idx};
            
            data.cueInfo(isnan([data.cueInfo.percent_respond])) = [];
            % data holds the bump data with bump mag as a variable
            
            mean_rt = []; std_rt = []; num_trials = []; freq_all = []; rt_data_all = []; freq = [];
            for d = 1:numel(data.cueInfo)
                mean_rt(d) = mean(data.cueInfo(d).rt);
                std_rt(d) = std(data.cueInfo(d).rt);
                num_trials(d) = numel(data.cueInfo(d).rt);
                freq(d) = frequencies_use(data.cueInfo(d).stimCode+1);
                if(use_std_err)
                    std_rt(d) = std_rt(d)/num_trials(d);
                end
                freq_all = [freq_all,freq(d)*ones(1,numel(data.cueInfo(d).rt))];
                rt_data_all = [rt_data_all,data.cueInfo(d).rt];
            end
            
            [fitObj,gof] = fit(freq_all',rt_data_all','a*exp(b*x)+c','startPoint',[0,0,0]);
            r2_all(counter,file_num) = gof.rsquare;
            
            x_data_fit = linspace(min(frequencies_use*0.9),max(frequencies_use)*1.1,100);
            y_data_fit = fitObj.a*exp(fitObj.b*x_data_fit)+fitObj.c;
            plot(x_data_fit,y_data_fit,'--','linewidth',2,'color',getColorFromList(1,color_idx(counter,file_num)));
            
            chan_all = (counter-1)+file_num + zeros(size(freq));
            monk_all = counter + zeros(size(freq));
            ANOVA_data.rt = [ANOVA_data.rt,mean_rt];
            ANOVA_data.freq = [ANOVA_data.freq, freq];
            ANOVA_data.monk = [ANOVA_data.monk, monk_all];
            ANOVA_data.chan = [ANOVA_data.chan, chan_all];
            
            % make plots
            % visual bar and horizontal line
            if(file_num==2)
                visual_bar = fill([0,max(frequencies_use)*2,max(frequencies_use)*2,0],...
                    mean_vis(1)+std_vis*[-1,-1,1,1],'k','linestyle','none');
                alpha(alpha_val);
                visual_line=plot([0,max(frequencies_use)*2],mean_vis+[0,0],'--','linewidth',2,'color','k');
                uistack(visual_line,'bottom');
                uistack(visual_bar,'bottom');
                            % bump bar and horz line
                if(~isempty(bump_data_rt))
                    bump_bar = fill([0,max(frequencies_use)*2,max(frequencies_use)*2,0],...
                        mean_bump(1)+std_bump*[-1,-1,1,1],'k','linestyle','none');
                    alpha(alpha_val);
                    bump_line=plot([0,max(frequencies_use)*2],mean_bump+[0,0],'-','linewidth',2,'color','k');
                    uistack(bump_line,'bottom');
                    uistack(bump_bar,'bottom');
                end
            end
            

            
            
            errorbar(freq+frequency_offset(file_num),mean_rt,std_rt,'.','color',getColorFromList(1,color_idx(counter,file_num)),...
                'markersize',28,'linewidth',2)
            hold on
            for d = 1:numel(data.cueInfo)
                dots_stim = scatter(freq(d)+frequency_offset(file_num)+zeros(1,numel(data.cueInfo(d).rt)),data.cueInfo(d).rt,...
                    28,'markeredgecolor','none','markerfacecolor',getColorFromList(1,color_idx(counter,file_num)));
                alpha(dots_stim, 0.5)
            end
        end
        % format figure
        xlim([0,max(frequencies_use)*1.1]);
        ylim([0,0.5]);
        formatForLee(gcf);
        xlabel('Frequency (Hz)');
        ylabel('RT (s)');
        
        set(gca,'fontsize',14)
        ax = gca;
        ax.XTick = [0:100:500];
        ax.XAxis.MinorTickValues = [0:50:500];

        counter = counter + 1;
    end
    
    
    
    ANOVA_data.chan = categorical(ANOVA_data.chan);
    ANOVA_data.monk = categorical(ANOVA_data.monk);
    
    ANOVA_data.chan = ANOVA_data.chan';
    ANOVA_data.rt = ANOVA_data.rt';
    ANOVA_data.monk = ANOVA_data.monk';
    ANOVA_data.freq = ANOVA_data.freq';
   %% do ANOVA
    lm_table = table(ANOVA_data.rt,ANOVA_data.chan,ANOVA_data.monk,ANOVA_data.freq,...
        'VariableNames',{'RT','chan','monk','freq'});
    
    lm = fitlm(lm_table,'RT~chan+monk+freq')