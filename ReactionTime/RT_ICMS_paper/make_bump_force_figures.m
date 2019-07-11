%% makes bump data figures
    monkey_names = {'Han','Duncan'};
    bump_data.folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Bump_data\';
    bump_data.alpha = 0.25;
    cd(bump_data.folderpath);
    use_std_err = 0; % otherwise will use standard deviations
    
    color_idx = [0,1];
    offset = 2*[-0.01,0.01;-0.05,0.05]/3;
%% plot bump staircase (2) for each monkey (2 figures total)
    
    counter = 1;
    r2_all = [];
    
    for monk = monkey_names
        file_list = dir([monk{1},'*']);
        f=figure();
        f.Name = [monk{1},'_bumpStaircase'];
        hold on
        for file_num = 1:numel(file_list)
            load(file_list(file_num).name);
            
            % sanitize data (remove stim trials)
            data.cueInfo([data.cueInfo.stimCode] ~= -1) = [];
            data.cueInfo(isEqual([data.cueInfo.bumpMag],1.1)) = [];
            
            % visual cue data
            vis_data_rt = data.cueInfo(isEqual([data.cueInfo.bumpMag],0)).rt;
            data.cueInfo(isEqual([data.cueInfo.bumpMag],0)) = [];
            mean_vis = mean(vis_data_rt);
            std_vis = std(vis_data_rt);
            num_vis_trials = numel(vis_data_rt);
            if(use_std_err)
                std_vis = std_vis./sqrt(num_vis_trials);
            end
            
            % data holds the bump data with bump mag as a variable
            bump_mags = [data.cueInfo.bumpMag];
            mean_rt = []; std_rt = []; num_trials = []; bump_mags_all = []; rt_data_all = [];
            for d = 1:numel(data.cueInfo)
                mean_rt(d) = mean(data.cueInfo(d).rt);
                std_rt(d) = std(data.cueInfo(d).rt);
                num_trials(d) = numel(data.cueInfo(d).rt);
                if(use_std_err)
                    std_rt(d) = std_rt(d)/sqrt(num_trials(d));
                end
                bump_mags_all = [bump_mags_all,bump_mags(d)*ones(1,numel(data.cueInfo(d).rt))];
                rt_data_all = [rt_data_all,data.cueInfo(d).rt];
            end
            
            [fitObj,gof] = fit(bump_mags_all',rt_data_all','a*exp(b*x)+c','startPoint',[0,0,0]);
            r2_all(counter,file_num) = gof.rsquare;
            
            x_data_fit = linspace(min(bump_mags*0.9),max(bump_mags)*1.1,100);
            y_data_fit = fitObj.a*exp(fitObj.b*x_data_fit)+fitObj.c;
            plot(x_data_fit,y_data_fit,'--','linewidth',2,'color',getColorFromList(1,color_idx(file_num)));
            
            
            
            % make plots
            % visual bar and horizontal line
            if(file_num==2)
                visual_bar = fill([0,max(bump_mags)*2,max(bump_mags)*2,0],...
                    mean_vis(1)+std_vis*[-1,-1,1,1],'k','linestyle','none');
                alpha(bump_data.alpha);
                visual_line=plot([0,max(bump_mags)*2],mean_vis+[0,0],'--','linewidth',2,'color','k');
                uistack(visual_line,'bottom');
                uistack(visual_bar,'bottom');
            end
            
            errorbar(bump_mags+offset(counter,file_num),mean_rt,std_rt,'.','color',getColorFromList(1,color_idx(file_num)),...
                'markersize',28,'linewidth',2)
            hold on
            for d = 1:numel(data.cueInfo)
                dots_bump = scatter(bump_mags(d)+offset(counter,file_num)+zeros(1,numel(data.cueInfo(d).rt)),data.cueInfo(d).rt,...
                    28,'markeredgecolor','none','markerfacecolor',getColorFromList(1,color_idx(file_num)))
                alpha(dots_bump,0.5)
            end
        end
        % format figure
        xlim([0,max(bump_mags)*1.1]);
        ylim([0,0.5]);
        formatForLee(gcf);
        xlabel('Perturbation force (N)');
        ylabel('RT (s)');
        
        set(gca,'fontsize',14)

        counter = counter + 1;
    end