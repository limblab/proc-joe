function [ outputData ] = getThreshold( outputData, inputData )
%
% preArtifact = 20;
% dataToTake = inputData.presample-preArtifact;
% 
% for art = 1:numel(outputData.artifactData)
%     outputData.artifactData(art).threshold = [];
%     data = outputData.artifactData(art).artifact;
%     meanData = repmat(mean(data(:,:,dataToTake),3),1,1,numel(dataToTake));
%     outputData.artifactData(art).threshold = rms(data(:,:,dataToTake)-meanData,3);
%     outputData.artifactData(art).artifact = data - repmat(mean(data(:,:,dataToTake),3),1,1,size(data,3));
% end

dataToTake = size(outputData.artifactData(1).artifact,3)-40:size(outputData.artifactData(1).artifact,3);

for art = 1:numel(outputData.artifactData)
    outputData.artifactData(art).threshold = [];
    data = outputData.artifactData(art).artifact;
    meanData = repmat(mean(data(:,:,dataToTake),3),1,1,numel(dataToTake));
    outputData.artifactData(art).threshold = rms(data(:,:,dataToTake)-meanData,3);
    outputData.artifactData(art).artifact = data - repmat(mean(data(:,:,dataToTake),3),1,1,size(data,3));
end

end

