function [ settlingTimeMetric ] = getSettlingTime( data, inputData )
% Words

settlingTimeMetric = zeros(size(data.artifactData.artifact,1),size(data.artifactData.artifact,2));
sampRate = 30000;
preSample = inputData.presample + 300e-6*sampRate;
% remove first 300 microseconds, then sum :D
for i = 1:size(data.artifactData.artifact,1)
    for j = 1:size(data.artifactData.artifact,2)
        settlingTimeMetric(i,j) = sum(abs(squeeze(data.artifactData.artifact(i,j,preSample:preSample+200))'))*1/sampRate;
    end
end

end

