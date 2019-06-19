%% makes 16 elec amp staircase data figures
    monkey_names = {'Han','Duncan'};
    map_files = {'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp',...
        'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\right S1 20180919\SN 6251-001804.cmp'};
    folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Multielec_staircase\';
    alpha_val = 0.25;
    cd(folderpath);
    use_std_err = 0; % otherwise will use standard deviations
    
    color_idx = [0,1];
    amp_offset = [-0.4,0.4];
    
%% plot stim staircase (train length, freq, amp) for each monkey (6 figures total)
    
    counter = 1;
    r2_all = [];
    
    for monk = monkey_names
    % train length data
        file_list = dir([monk{1},'*EXAMPLE*data*']);
        
        f=figure();
        f.Name = [monk{1},'_16elecs_staircase'];
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
                stim_amps_all = [stim_amps_all,STIM_PARAMS(d)*ones(1,numel(data.cueInfo(d).rt))];
                rt_data_all = [rt_data_all,data.cueInfo(d).rt];
            end
            
%             if(strcmpi(monk{1},'Han') == 1)
%                 [fitObj,gof] = fit(stim_amps_all',rt_data_all','a*exp(b*x)+c','startPoint',[0,0,0]);
%                 r2_all(counter,file_num) = gof.rsquare;
% 
%                 x_data_fit = linspace(min(STIM_PARAMS*0.9),max(STIM_PARAMS)*1.1,100);
%                 y_data_fit = fitObj.a*exp(fitObj.b*x_data_fit)+fitObj.c;
%                 plot(x_data_fit,y_data_fit,'--','linewidth',2,'color',getColorFromList(1,color_idx(file_num)));
%             end
 
            % make plots
            % visual bar and horizontal line
            if(~isempty(vis_data_rt))
                visual_bar = fill([0,max(STIM_PARAMS)*2,max(STIM_PARAMS)*2,0],...
                    mean_vis(1)+std_vis*[-1,-1,1,1],'k','linestyle','none');
                alpha(visual_bar,alpha_val);
                visual_line=plot([0,max(STIM_PARAMS)*2],mean_vis+[0,0],'--','linewidth',2,'color','k');
                uistack(visual_line,'bottom');
                uistack(visual_bar,'bottom');
            end
            
            % bump bar and horz line
            if(~isempty(bump_data_rt))
                bump_bar = fill([0,max(STIM_PARAMS)*2,max(STIM_PARAMS)*2,0],...
                    mean_bump(1)+std_bump*[-1,-1,1,1],'k','linestyle','none');
                alpha(bump_bar,alpha_val);
                bump_line=plot([0,max(STIM_PARAMS)*2],mean_bump+[0,0],'-','linewidth',2,'color','k');
                uistack(bump_line,'bottom');
                uistack(bump_bar,'bottom');
            end
            
            errorbar(STIM_PARAMS+amp_offset(file_num),mean_rt,std_rt,'.','color',getColorFromList(1,color_idx(file_num)),...
                'markersize',28,'linewidth',2)
            hold on
            for d = 1:numel(data.cueInfo)
                dots_stim = scatter(STIM_PARAMS(d)+amp_offset(file_num)+zeros(1,numel(data.cueInfo(d).rt)),data.cueInfo(d).rt,...
                    28,'markeredgecolor','none','markerfacecolor',getColorFromList(1,color_idx(file_num)))
                alpha(dots_stim,0.5)
            end
        end
        % format figure
        xlim([0,max(STIM_PARAMS)*1.1]);
        ylim([0,0.5]);
        formatForLee(gcf);
        xlabel('Amplitude (\muA)');
        ylabel('RT (s)');
        
        set(gca,'fontsize',14)
        ax = gca;
        ax.XTick = [0,STIM_PARAMS];
        ax.XAxis.MinorTickValues = [0:5:100];
        

        counter = counter + 1;
    end
    
    
%% make array pics with patterns colored in
    counter = 1;
    for monk = monkey_names
        f=figure();
        f.Name = [monk{1},'_16elec_stimStaircase_arrayMap'];
        % get pattern data and map file data
        file_list = dir([monk{1},'*EXAMPLE*elecList*']);
        load(file_list.name);
        map_data = loadMapFile(map_files{counter});
        
        for p = 1:2
            switch p
                case 1
                    chan_list = pattern1;
                case 2
                    chan_list = pattern2;
            end
                
            hold on
            for elec = 1:numel(chan_list)
                map_idx = find(chan_list(elec) == map_data.chan);
                rectangle('Position',[map_data.row(map_idx),11-map_data.col(map_idx),1,1],...
                    'FaceColor',getColorFromList(1,color_idx(p)), 'EdgeColor','none');
                hold on
            end
        end
        plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5)

        set(gca,'visible','off')
        xlim([1,11])
        ylim([1,11])
        axis square
        
        counter = counter + 1;
    end
