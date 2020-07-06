function [output_data] = stratifySamplePDs(data, n_groups)
    % group data groups of size n_groups based on proximity. Minimize MSE within each group.
    % Sample an electrode from each group to generate n_groups. Excess data
    % (groups that are not full) will not be used. 

    % useful variables
    n_data = numel(data);
    groups = nan(floor(n_data/n_groups),n_groups);
    n_excess = n_data - floor(n_data/n_groups)*n_groups;
    group_size = n_groups;
    
    %% group data into groups of size group_size
    % order data
    [data,sort_idx] = sort(data);
    
    % PD data is circular, so all we need to do is shift the starting
    % point, group and compute an MSE. Then pick best MSE
    best_MSE = 1000000; % set as very large
    best_group_idx = -1; best_excess_idx = -1; % will be overriden
    % get list of idx to test removing
    remove_data = combnk((1:1:numel(data)),n_excess);
    % test all remove, and start positions
    for remove_idx = 1:size(remove_data,1)
        % remove excess data
        data_trim = data; data_trim(remove_data(remove_idx,:)) = [];
        % setup group_idx
        group_idx = nan(numel(data_trim),1);
        group_idx(1:numel(data_trim)) = ceil([1/group_size:1/group_size:(numel(data_trim))/group_size]);
        for start_idx = 1:group_size % consider each possible starting spot
            % compute MSE assigning data to group_idx
            curr_MSE = 0;
            for i_group = 1:size(groups,1) % for each group
                mean_group = circular_mean(data_trim(group_idx == i_group));
                group_MSE = circular_mean(circ_dist(data_trim(group_idx == i_group),mean_group).^2);
                curr_MSE = curr_MSE + group_MSE;
            end

            % determine if better than previous groups
            if(curr_MSE < best_MSE)
                best_MSE = curr_MSE;
                best_group_idx = group_idx;
                best_excess_idx = remove_data(remove_idx,:);
            end

            % shift group idx by one
            group_idx = circshift(group_idx,1);
        end
    end
    
    % add nan's at remove_idx in group_idx
    for i_remove = numel(best_excess_idx):-1:1
        if(best_excess_idx(i_remove) >= numel(group_idx))
            group_idx = [group_idx;nan];
        elseif(best_excess_idx(i_remove) ==1)
            group_idx = [nan;group_idx];
        else
            group_idx = [group_idx(1:best_excess_idx(i_remove)-1); nan; group_idx(best_excess_idx(i_remove):end)];
        end
    end
    
    % sample once from each group to get groups
    group_id = [1:1:group_size];
    for i_group = 1:size(groups,1) % for each group
        % shuffle group_id
        group_id = group_id(randperm(numel(group_id)));
        
        % place data point into corresponding group in groups
        group_data = sort_idx(best_group_idx == i_group);
        groups(i_group,:) = group_data(group_id); % place into groups
    end
    
    
    % package outputs
    output_data.groups = groups;
    output_data.MSE = best_MSE;
    output_data.group_idx = group_idx;
    output_data.sort_idx = sort_idx;
end


function r =  circ_dist(x,y)
%
% r = circ_dist(alpha, beta)
%   Pairwise difference x_i-y_i around the circle computed efficiently.
%
%   Input:
%     alpha      sample of linear random variable
%     beta       sample of linear random variable or one single angle
%
%   Output:
%     r       matrix with differences
%
% References:
%     Biostatistical Analysis, J. H. Zar, p. 651
%
% PHB 3/19/2009
%
% Circular Statistics Toolbox for Matlab
% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
    if size(x,1)~=size(y,1) && size(x,2)~=size(y,2) && length(y)~=1
      error('Input dimensions do not match.')
    end
    r = angle(exp(1i*x)./exp(1i*y));
end