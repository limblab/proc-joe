%% make matrix of response (bins) by num neurons x num channels, for each condition 

    pca_data = [];
    pca_data.data = cell(numel(arrayData{1}.STIM_PARAMETERS),1);
    for cond = 1:numel(arrayData{1}.STIM_PARAMETERS)
        for arr_idx = 1:numel(arrayData)
            for chan = 1:size(arrayData{arr_idx}.bC,1)
                pca_data.data{cond} = [pca_data.data{cond}, arrayData{arr_idx}.bC{chan,cond}'];
            end
        end
    end


%%  do pca on each condition

    for cond = 4%:numel(arrayData{1}.STIM_PARAMETERS)
        [coeff,score,latent,tsquared,explained,mu] = pca(pca_data.data{cond}(80:180,:));
        
        
    end
