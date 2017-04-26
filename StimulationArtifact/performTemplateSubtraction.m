function [ outputData ] = performTemplateSubtraction( outputData )
% performs template subtraction, does this on cathodal and anodal pulses
% separately (cathodal = odd, anodal = even)

for art = 1:numel(outputData.artifactData)
    data = outputData.artifactData(art).artifact;
    templateCathodal = mean(data(:,1:2:end,:),2);
    templateAnodal = mean(data(:,2:2:end,:),2); % mean across the 100 electrodes, results in a 96x1xnumPoints template
    
    data(:,1:2:end,:) = data(:,1:2:end,:)-repmat(templateCathodal,1,size(data,2)/2,1);
    data(:,2:2:end,:) = data(:,2:2:end,:)-repmat(templateAnodal,1,size(data,2)/2,1);
    outputData.artifactData(art).artifact = data; 
end

end

