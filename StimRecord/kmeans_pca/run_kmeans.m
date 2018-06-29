function [ outputData ] = run_kmeans( data,opts )
    % data is an N x M matrix where N = # observations (data points), M = #
    % variables for that data

    %% configure opts
    opts = configureOpts(opts);
    if(isempty(opts.INDEX))
        opts.INDEX = [1,size(data,2)];
    end
    
    %% initialize everything
    cluster = zeros(opts.MAX_ITER,size(data,1));
    cluster_center = datasample(data,opts.NUM_CLUSTERS,'replace',false);
    
    %% subtract baseline from data
    if(opts.SUBTRACT_BASELINE)
        data = data - mean(data,2);
    end
    
    %% run k-means
    for iter = 1:opts.MAX_ITER
        % assign data to a cluster
        dist = computeDistance(data,cluster_center,opts);
        [~,cluster(iter,:)] = min(dist,[],2);
        
        % move cluster centers
        for c = 1:numel(opts.NUM_CLUSTERS)
            mask = cluster(iter,:) == c;
            if(sum(mask)~=0)
                cluster_center(c,:) = mean(data(mask,:));
            end
        end
    end
    
    %% output things
    outputData.cluster = cluster;
    outputData.cluster_center = cluster_center;
    
end


function [dist] = computeDistance(data,cluster_center,opts)


    % reformat data and cluster_center
    data  = data(:,opts.INDEX(1):opts.INDEX(2));
    cluster_center = cluster_center(:,opts.INDEX(1):opts.INDEX(2));
    
    
    data = permute(data,[1,3,2]);
    data = repmat(data,1,opts.NUM_CLUSTERS,1);
    
    cluster_center = permute(cluster_center,[3,1,2]);
    cluster_center = repmat(cluster_center,size(data,1),1,1);
    % compute dist
    if(strcmpi(opts.DISTANCE_METHOD,'Euclidean')==1)
        dist = sum(((data-cluster_center).^2),3);
    else % default is a cosine similarity
        dist = 1-dot(data,cluster_center,3)./sqrt(dot(data,data,3).*dot(cluster_center,cluster_center,3));
    end
end

function [opts] = configureOpts(optsInput)

    opts.NUM_CLUSTERS = 3;
    opts.MAX_ITER = 100;
    opts.DISTANCE_METHOD = 'cosine'; % default case
    opts.INDEX = [];
    opts.SUBTRACT_BASELINE = 1;
    
    %% check if in optsSave and optsSaveInput, overwrite if so
    try
        inputFieldnames = fieldnames(optsInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(opts,inputFieldnames{fn}))
               opts.(inputFieldnames{fn}) = optsInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
    
end