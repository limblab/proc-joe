function [outputData] = performPCAStep(outputData,coeffsRemove)
% this performs PCA to remove artifacts using the first component only

for art = 1:numel(outputData.artifactData)
    for stim = 1:size(outputData.artifactData(art).artifact,2)
%     for elec = 1:size(outputData.artifactData(art).artifact,1)
        % run PCA
        [coeff,score,latent,tsquared,explained,mu] = pca(squeeze(outputData.artifactData(art).artifact(:,stim,:)),'Centered',false);
%         [coeff,score,latent,tsquared,explained,mu] = pca(squeeze(outputData.artifactData(art).artifact(elec,:,:)),'Centered',false);
        outputData.coeff = coeff;
        outputData.score = score;
        coeff(:,1:coeffsRemove) = 0;
        artifactPCA = score*coeff';
        outputData.artifactData(art).artifact(:,stim,:) = artifactPCA;
%         outputData.artifactData(art).artifact(elec,:,:) = artifactPCA;
    end
end

end