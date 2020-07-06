function [space_constant, fit_data] = getSpaceConstant(resp, dist, dist_bin_edges)
    % define outputs
    space_constant = [];
    fit_data = [];

    
    % bin data
    [bin_dist,~,bin_idx] = histcounts(dist,dist_bin_edges);
    bin_centers = dist_bin_edges(1:end-1) + mode(diff(dist_bin_edges))/2;
    
    bin_resp = []; bin_dist = [];
    for i_bin = 1:numel(bin_centers)
        if(sum(bin_idx == i_bin) > 0)
            bin_resp(end+1,1) = mean(resp(bin_idx == i_bin));
            bin_dist(end+1,1) = bin_centers(i_bin);
        end
    end
            
    % fit data
    [fit_data.f,fit_data.gof] = fit(bin_dist,bin_resp,...
        'a*exp(-x/b)+c','Lower',[0,0,-0.4],'Upper',[5,10000,0.4],'StartPoint',[0.03,-0.005,0]);

    % get space constant
    if(fit_data.gof.rsquare > 0.1)
        space_constant = fit_data.f.b;
    else
        space_constant = nan;
    end
    
end