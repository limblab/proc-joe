function [outputData] = performCurveFittingStep(outputData)
% this function performs curve fitting on the data provided

% curve fitting on all channels
startIdx = 13; % be more sophisticated with this
endIdx = 30; % be more sophisticated with this

for art = 1:numel(outputData.artifactData)
    for ch = 1:size(outputData.artifactData(art).artifact,1)
        data = squeeze(outputData.artifactData(art).artifact(ch,1:2:end,startIdx:endIdx))'; % 50 x 280ish matrix
        sizeData = size(data);
        data = reshape(data,numel(data),1);
        x = repmat((startIdx:1:endIdx)',length(data)/(endIdx-startIdx+1),1);
        expEqn = 'poly4';
        % fit data with expEqun

        f=fit(x,data,expEqn);
        fitData = data - feval(f,x);
        outputData.artifactData(art).artifact(ch,1:2:end,startIdx:endIdx) = reshape(fitData,sizeData(1),sizeData(2))';
        
        data = squeeze(outputData.artifactData(art).artifact(ch,2:2:end,startIdx:endIdx))'; % 50 x 280ish matrix
        sizeData = size(data);
        data = reshape(data,numel(data),1);
        x = repmat((startIdx:1:endIdx)',length(data)/(endIdx-startIdx+1),1);
        expEqn = 'poly4';
        % fit data with expEqun

        f=fit(x,data,expEqn);
        fitData = data - feval(f,x);
        outputData.artifactData(art).artifact(ch,2:2:end,startIdx:endIdx) = reshape(fitData,sizeData(1),sizeData(2))';
    end
end

end