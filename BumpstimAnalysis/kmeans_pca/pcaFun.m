%% make a matrix of data: N x M where N = number of units/stim observations
% and M = data points

mat = [];
for arr = 1:numel(arrayData)
    for ampIdx = 1:numel(arrayData{arr}.spikeTrialTimes)
        mat(end+1,:) = arrayData{arr}.bC{ampIdx};
    end
end
xData = (arrayData{1}.bE{1}(1:end-1) + mode(diff(arrayData{1}.bE{1}))/2)';

% run pca

[coeff,score,latent,tsquared,explained,mu] = pca(mat','NumComponents',4);

% plot components
figure();
plot(xData,score,'linewidth',1.5)


%% run k-means to categorize responses

opts.NUM_CLUSTERS = 3;
opts.MAX_ITER = 500;
opts.INDEX = [40,50];
opts.DISTANCE_METHOD = 'Euclidean';
[clusterData] = run_kmeans(mat,opts);

% plot each clusters' waveforms
figure();
for c = 1:opts.NUM_CLUSTERS
    subplot(opts.NUM_CLUSTERS,1,c);
    mask = clusterData.cluster(end,:) == c;
    if(sum(mask)~=0)
        plot(xData',mat(mask,:))
    end
end

% plot each cluster's center
colors = {'r',[0,0.5,0],'b'};
figure();
for c = 1:opts.NUM_CLUSTERS
    plot(xData,clusterData.cluster_center(c,:),'color',colors{c})
    hold on
end
% add cluster info to arrayData
counter = 1;
for arr = 1:numel(arrayData)
    for ampIdx = 1:numel(arrayData{arr}.spikeTrialTimes)
        arrayData{arr}.cluster(ampIdx,1) = clusterData.cluster(end,counter);
        counter = counter + 1;
    end
end
%% plot array with cluster information?

for ampIdx = 1:numel(arrayData{arr}.spikeTrialTimes)
    figure();
    
    plot([1,11,11,1,1],[1,1,11,11,1],'-k','linewidth',1.5)
    ax = gca;
    ax.YTickLabel = {};
    ax.XTickLabel = {};
    for arr = 1:numel(arrayData)
        rectangle('Position',[arrayData{arr}.COL,arrayData{arr}.ROW,1,1],'EdgeColor','k',...
                'FaceColor',colors{arrayData{arr}.cluster(ampIdx)},'linewidth',0.1);
    end
    set(gca,'Visible','off');
    axis square
end
