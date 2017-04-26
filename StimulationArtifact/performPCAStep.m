function [outputData] = performPCAStep(outputData)
% this performs PCA to remove artifacts using the first component only

for art = 1:numel(outputData.artifactData)
    for stim = 1:size(outputData.artifactData(art).artifact,2)
        % run PCA
        [coeff,score,latent,tsquared,explained,mu] = pca(squeeze(outputData.artifactData(art).artifact(:,stim,:)),'Centered',false);
        coeff(:,1:1) = 0;
        artifactPCA = score*coeff';
        outputData.artifactData(art).artifact(:,stim,:) = artifactPCA;
    end
end

end