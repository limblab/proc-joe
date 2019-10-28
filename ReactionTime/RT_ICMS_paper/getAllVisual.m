%% all visual data

    folderpaths = {'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Bump_data\';...
        'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\SingleElecs\';...
        'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\SingleElecs_staircase\';...
        'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Multielec_staircase\';...
        'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\Distance_data\';...
        'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\SingleElecs_200uA\';...
        'C:\Users\Joseph\Desktop\Lab\Data\ReactionTime\RandomElecs\'};
    
    vis_mean_aggregate = [];
    vis_std_aggregate = [];
    vis_num_trials_aggregate = [];
    vis_monkey_name_aggregate = {};
    vis_date_aggregate = {};
    
    
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

                data.cueInfo([data.cueInfo.bumpMag]>0) = [];

                vis_idx = [];
                vis_idx = find([data.cueInfo.bumpMag] == 0);


                % data should just have vis data
                if(~isempty(vis_idx) && ~isempty(data.cueInfo(vis_idx).rt))
                    vis_mean_aggregate(end+1,1) = mean(data.cueInfo(vis_idx).rt);
                    vis_std_aggregate(end+1,1) = std(data.cueInfo(vis_idx).rt);
                    vis_num_trials_aggregate(end+1,1) = numel(data.cueInfo(vis_idx).rt);
                    vis_monkey_name_aggregate{end+1} = monkey_name;
                else
                    warning(['No vis data for file ',file_list(file).name]);
                end
            catch
                warning(['No vis data for file ',file_list(file).name]);
            end
        end
    end
    
    
        
    is_han = strcmpi(vis_monkey_name_aggregate,'Han')==1;
    han_vis_pool_std = pooled_std(vis_std_aggregate(is_han),vis_num_trials_aggregate(is_han));
    han_vis_mean = mean(vis_mean_aggregate(is_han));
    duncan_vis_pool_std = pooled_std(vis_std_aggregate(~is_han),vis_num_trials_aggregate(~is_han));
    duncan_vis_mean = mean(vis_mean_aggregate(~is_han));
