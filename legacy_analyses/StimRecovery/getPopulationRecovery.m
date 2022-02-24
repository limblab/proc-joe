function [recov_data] = getPopulationRecovery(neural_act, bin_edges)

    % finds time to recovery after stim onset (0ms in bin_edges)
    bin_centers = bin_edges(1:end-1) + mode(diff(bin_edges))/2;

    
    % use mean neural act across stims to compute this
    mean_neural_act = squeeze(mean(neural_act,2));
        
    % compute recovery for entire space
    [all_dim_recov,all_dim_fit,dist_to_baseline] = getRecoveryHelper(mean_neural_act, bin_centers);
    
    % compute recovery for each dimension
%     for i_dim = 1:size(mean_neural_act,1)
%         single_dim_recov(i_dim) = getRecoveryHelper(mean_neural_act(i_dim,:), bin_centers);
%     end


    recov_data.all_dim = all_dim_recov;
    recov_data.dist_to_baseline = dist_to_baseline;
    recov_data.all_dim_fit = all_dim_fit;
    
%     recov_data.single_dim = single_dim_recov;
end



function [recov_time,fit_data,dist_to_baseline] = getRecoveryHelper(mean_neural_act, bin_centers)

    % get baseline vals and post stim vals
    baseline_idx = find(bin_centers < 0);
    baseline_idx(end) = [];
    
    mean_neural_act = mean_neural_act - mean(mean_neural_act(:,baseline_idx),2);
    mean_baseline = mean(mean_neural_act(:,baseline_idx),2);
    
    % get distance to baseline vals across all dims provided
    dist_to_baseline = sqrt(sum((mean_neural_act - mean_baseline).^2,1));
    thresh = mean(dist_to_baseline(baseline_idx)) + std(dist_to_baseline(baseline_idx));    
    
    % define time recovered from stim onset (0ms) based on distance to
    % baseline vals across all dims
    post_stim_idx = find(bin_centers > 0);
    post_stim_offset_idx = find(bin_centers > 100,1,'first');
    dist_post_stim = dist_to_baseline(post_stim_idx);
    
    idx_recov = find(dist_post_stim < thresh,1,'first');
    
    if(isempty(idx_recov))
        idx_recov = nan;
        end
    recov_time = idx_recov*mode(diff(bin_centers)) - mode(diff(bin_centers))/2;
    
    
    % fit exponential decay
    y_data = dist_to_baseline(post_stim_offset_idx:1:end) - mean(dist_to_baseline(1:post_stim_idx-1));
    x_data = ((1:1:numel(y_data))-1)*mode(diff(bin_centers));
    
    fit_data = fit(x_data',y_data','a*exp(-b*x)','Lower',[0,0],'Upper',[50,1000]);
    


end