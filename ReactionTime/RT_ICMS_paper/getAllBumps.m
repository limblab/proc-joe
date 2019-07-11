%% all_bump data

    folderpaths = {'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Bump_data\';...
        'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\SingleElecs\';...
        'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\SingleElecs_staircase\';...
        'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Multielec_staircase\';...
        'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Distance_data\';...
        'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\SingleElecs_200uA\';...
        'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\RandomElecs\'};
    
    bump_mean_aggregate = [];
    bump_std_aggregate = [];
    bump_num_trials_aggregate = [];
    bump_monkey_name_aggregate = {};
    p_vals_bump_vis = [];
    
    for fold = folderpaths'
        disp(fold)
        cd(fold{1});
        file_list = dir('*data*');
        
        for file = 1:numel(file_list)
            
            try
                load([fold{1},file_list(file).name]);

                % get monkey name
                underscore_idx = strfind(file_list(file).name,'_');
                monkey_name = file_list(file).name(1:underscore_idx(1)-1);
                

                 % sanitize data (remove stim trials and visual trials)
                data.cueInfo([data.cueInfo.stimCode] ~= -1) = [];

                data.cueInfo(isEqual([data.cueInfo.bumpMag],1.1)) = [];
                
                vis_data = data.cueInfo(isEqual([data.cueInfo.bumpMag],0));
                data.cueInfo(isEqual([data.cueInfo.bumpMag],0)) = [];
        
                bump_idx = [];
                if(strcmpi(monkey_name,'Han'))
                    bump_idx = find([data.cueInfo.bumpMag] == 1);
                elseif(strcmpi(monkey_name,'Duncan'))
                    bump_idx = find([data.cueInfo.bumpMag] == 3);
                end

                if(~isempty(bump_idx) && ~isempty(vis_data.rt))
                    rt_vis = vis_data.rt;
                    rt_bump = data.cueInfo(bump_idx).rt;
                    tail = 'left';
                    [h,p,ci,stats] = ttest2(rt_bump,rt_vis,0.95,tail,'unequal');
                    p_vals_bump_vis(end+1,1) = p;
                    if(numel(p_vals_bump_vis) == 7)
                        disp('remember to remove the outlier');
                    end
                end
                % data should just have bump data
                if(~isempty(bump_idx))
                    bump_mean_aggregate(end+1,1) = mean(data.cueInfo(bump_idx).rt);
                    bump_std_aggregate(end+1,1) = std(data.cueInfo(bump_idx).rt);
                    bump_num_trials_aggregate(end+1,1) = numel(data.cueInfo(bump_idx).rt);
                    bump_monkey_name_aggregate{end+1} = monkey_name;
                    bump_file_name_aggregate{end+1} = file_list(file).name;
                end
            catch
                warning(['No bump data for file ',file_list(file).name]);
            end
        end
    end
    
    is_han = strcmpi(bump_monkey_name_aggregate,'Han')==1;
    han_bump_pool_std = pooled_std(bump_std_aggregate(is_han),bump_num_trials_aggregate(is_han));
    han_bump_mean = mean(bump_mean_aggregate(is_han));
    duncan_bump_pool_std = pooled_std(bump_std_aggregate(~is_han),bump_num_trials_aggregate(~is_han));
    duncan_bump_mean = mean(bump_mean_aggregate(~is_han));
    
