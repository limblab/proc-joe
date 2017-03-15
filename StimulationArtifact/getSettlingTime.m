function [ settlingTimeMetric ] = getSettlingTime( data, inputData )
% Words

settlingTimeMetric = zeros(numel(data.artifactData),size(data.artifactData(1).artifact,1),size(data.artifactData(1).artifact,2));
sampRate = 30000;
preSample = inputData.presample + 300e-6*sampRate;
% remove first 300 microseconds, then sum :D
for a = 1:numel(data.artifactData)
    for i = 1:size(data.artifactData(a).artifact,1)
        for j = 1:size(data.artifactData(a).artifact,2)
            settlingTimeMetric(a,i,j) = sum(abs(squeeze(data.artifactData(a).artifact(i,j,preSample:preSample+inputData.windowSize-1-300e-6*sampRate))'))*1/sampRate;
        end
    end
end

% compress settlingTimeMetric across first dim
settlingTimeMetric = squeeze(mean(settlingTimeMetric,1));
end

