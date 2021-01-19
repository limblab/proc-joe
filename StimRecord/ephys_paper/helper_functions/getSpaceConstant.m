function [space_constant, rsquare] = getSpaceConstant(resp, dist, dist_bin_edges)
    % define outputs
    space_constant = [];
    rsquare = [];
    
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
    try
        [fit_data.f,fit_data.gof] = fit(bin_dist,bin_resp,...
            'a*exp(-x/b)+c','Lower',[0,0,-0.4],'Upper',[1,7500,1],'StartPoint',[0.5,1000,0]);

        % get space constant
        space_constant = fit_data.f.b;
        rsquare = fit_data.gof.rsquare;
    catch
        space_constant = nan;
        rsquare = nan;
    end
    
end