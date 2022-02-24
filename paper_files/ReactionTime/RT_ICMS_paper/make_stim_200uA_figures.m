%% makes stim 200uA staircase data figures
    monkey_names = {'Han','Duncan'};
    folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\SingleElecs_200uA\';
    alpha_val = 0.25;
    cd(folderpath);
    use_std_err = 0; % otherwise will use standard deviations
    
    color_idx = [4,5]; %chan 19 first, then chan 42 for Duncan
    amp_offset = [-1.5,1.5];
    stim_amps = [40:40:200];
    
%% plot stim staircase for each monkey up to 200uA (2 figures)
    
    counter = 1;
    r2_all = [];
    
    for monk = monkey_names(1)
        file_list = dir([monk{1},'*EXAMPLE*']);
        
        f=figure();
        f.Name = [monk{1},'_singleElecs_200uA'];
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
            
            % data holds the bump data with bump mag as a variable
            mean_rt = []; std_rt = []; num_trials = []; stim_amps_all = []; rt_data_all = [];
            
            for d = 1:numel(data.cueInfo)
                mean_rt(d) = mean(data.cueInfo(d).rt);
                std_rt(d) = std(data.cueInfo(d).rt);
                num_trials(d) = numel(data.cueInfo(d).rt);
                if(use_std_err)
                    std_rt(d) = std_rt(d)/sqrt(num_trials(d));
                end
                stim_amps_all = [stim_amps_all,stim_amps(d)*ones(1,numel(data.cueInfo(d).rt))];
                rt_data_all = [rt_data_all,data.cueInfo(d).rt];
            end
            
            [fitObj,gof] = fit(stim_amps_all',rt_data_all','a*exp(b*x)+c','startPoint',[0,0,0]);
            r2_all(counter,file_num) = gof.rsquare;
            
            x_data_fit = linspace(min(stim_amps*0.9),max(stim_amps)*1.1,100);
            y_data_fit = fitObj.a*exp(fitObj.b*x_data_fit)+fitObj.c;
            plot(x_data_fit,y_data_fit,'--','linewidth',2,'color',getColorFromList(1,color_idx(file_num)));
            
            num_trials_all = [num_trials_all,num_trials,vis_num_trials,bump_num_trials];
 
            % make plots
            % visual bar and horizontal line
            if(~isempty(vis_data_rt))
                visual_bar = fill([0,max(stim_amps)*2,max(stim_amps)*2,0],...
                    mean_vis(1)+std_vis*[-1,-1,1,1],'k','linestyle','none');
                alpha(alpha_val);
                visual_line=plot([0,max(stim_amps)*2],mean_vis+[0,0],'--','linewidth',2,'color','k');
                uistack(visual_line,'bottom');
                uistack(visual_bar,'bottom');
            end
            
            % bump bar and horz line
            if(~isempty(bump_data_rt))
                bump_bar = fill([0,max(stim_amps)*2,max(stim_amps)*2,0],...
                    mean_bump(1)+std_bump*[-1,-1,1,1],'k','linestyle','none');
                alpha(alpha_val);
                bump_line=plot([0,max(stim_amps)*2],mean_bump+[0,0],'-','linewidth',2,'color','k');
                uistack(bump_line,'bottom');
                uistack(bump_bar,'bottom');
            end
            
            errorbar(stim_amps+amp_offset(file_num),mean_rt,std_rt,'.','color',getColorFromList(1,color_idx(file_num)),...
                'markersize',28,'linewidth',2)
            hold on
            for d = 1:numel(data.cueInfo)
                dots_stim = scatter(stim_amps(d)+amp_offset(file_num)+zeros(1,numel(data.cueInfo(d).rt)),data.cueInfo(d).rt,...
                    28,'markeredgecolor','none','markerfacecolor',getColorFromList(1,color_idx(file_num)))
                alpha(dots_stim,0.5)
            end
        end
        % format figure
        xlim([0,max(stim_amps)+20]);
        ylim([0,0.5]);
        formatForLee(gcf);
        xlabel('Amplitude (\muA)');
        ylabel('RT (s)');
        
        set(gca,'fontsize',14)
        ax = gca;
        ax.XTick = [0,stim_amps];
        ax.XAxis.MinorTickValues = [0:20:200];
        

        counter = counter + 1;
    end