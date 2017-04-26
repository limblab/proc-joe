function [ settlingTimeMetric ] = getSettlingTime( data, inputData )
% Words

settlingTimeMetric = zeros(numel(data.artifactData),size(data.artifactData(1).artifact,1),size(data.artifactData(1).artifact,2));
sampRate = 30000;
preSample = inputData.presample + 300e-6*30000; % remove 300 microseconds
windowSize = inputData.windowSize;
% go backwards in time, find first time n out of m data points are out of
% desired range: mean + 2*std?
% for a = 1:numel(data.artifactData)
%     for i = 1:size(data.artifactData(a).artifact,1)
%         for j = 1:size(data.artifactData(a).artifact,2)
%             if(sum(i-inputData.badChList==0)==0) % not in bad ch list (inputData.badChList)
%                 dataToLookAt = squeeze(data.artifactData(a).artifact(i,j,preSample:end));
%                 settlingTimeMetric(a,i,j) = computeSettlingTime(dataToLookAt,sampRate,windowSize); 
%             end
%         end   
%     end  
% end

% summation idea below
% remove first 300 microseconds, then sum :D
for a = 1:numel(data.artifactData)
    for i = 1:size(data.artifactData(a).artifact,1)
        for j = 1:size(data.artifactData(a).artifact,2)
            dataToAnalyze = squeeze(data.artifactData(a).artifact(i,j,preSample:end));
            meanDataToAnalyze = mean(dataToAnalyze(numel(dataToAnalyze)/2:end));
            settlingTimeMetric(a,i,j) = sum(abs(dataToAnalyze-meanDataToAnalyze))*1/sampRate;
        end
    end
end

% compress settlingTimeMetric across first dim
settlingTimeMetric = squeeze(mean(settlingTimeMetric,1));
end

