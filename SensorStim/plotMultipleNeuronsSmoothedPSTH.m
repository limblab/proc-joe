% plots a PSTH a single for every muscle that was stimulated
% this script assumes variable 'common' is in workspace

gaussianSmooth_std = 0.05*2; % 50 ms bin widths generally used
binSize = 0.05;
nn=1; % neuron number in common

% % perform gaussian smoothing?
kernel_width = ceil(3*gaussianSmooth_std/binSize);
kernel = normpdf(-kernel_width*binSize: ...
    binSize:kernel_width*binSize,...
    0, gaussianSmooth_std); 
normalizer = conv(kernel,ones(1,length(binEdges)-1));
for i = 1:size(neuronCounts,1)
    smoothed_fr_inter = conv(kernel,neuronCounts(i,:))./normalizer;
    neuronCounts(i,:) = smoothed_fr_inter(kernel_width+1:end-kernel_width);
end

% plot neuron for all muscles
figure();
hold on
for i = 1:length(common.muscleIdx)
    yData = common.(genvarname(common.muscleIdx{i,1})).counts{nn};
    xData = common.(genvarname(common.muscleIdx{i,1})).edges{nn};
    
    kernel_width = ceil(3*gaussianSmooth_std/binSize);
    kernel = normpdf(-kernel_width*binSize: ...
    binSize:kernel_width*binSize,...
    0, gaussianSmooth_std); 
    normalizer = conv(kernel,ones(1,length(xData)-1));
    smoothed_fr_inter = conv(kernel,yData)./normalizer;
    yData = smoothed_fr_inter(kernel_width+1:end-kernel_width);

    yData = yData-sum(yData)/length(yData);
    yData = yData/max(yData);
    plot(xData(1:end-1)+mode(diff(xData))/2,yData);
end
