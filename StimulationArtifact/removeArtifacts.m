function [outputDataFiltered] = removeArtifacts(outputData,inputData,removalSteps,folderpath)
% this function removes artifacts according to removalSteps.
% outputDataFiltered{idx} corresponds to removalSteps(2,:)
% removal steps is encoded as follows: 
%   removalSteps{~,1} = order of steps -- 0 = filter, 1 = curve fit, 2 =
%   PCA. so '012' would mean filter, then curve fit, then do PCA. '01'
%   would be filter, curve fit
%   removalSteps{~,2} is information related to the filtering step. 
%   Filter information: low pass  = '*Low#', high pass = '*High#', order =
%   '*Order#'... example 'Low1000_High10_Order6' would make a 6th order
%   filter with a passband between 10Hz and 1000Hz. If no information is
%   provided, -1 is the default corresponding to no low pass/no high pass. 
%   removalSteps{~,3} is information related to the filtering step.  Curve 
%   fit information: 
%   removalSteps{~,4} is information related to the pca step PCA information: 
warning('off')

outputDataFiltered = {};
outputData = getThreshold(outputData,inputData);

for removalStepIdx = 1:size(removalSteps,1)
    % init outputDataTemp
    outputDataTemp = outputData;
    % get order of steps
    stepOrder = removalSteps{removalStepIdx,:};
   
    % move through step order to remove artifacts
    for step = 1:length(stepOrder)
        if(stepOrder(step) == '0') % filtering step
            lowPass = -1;
            highPass = -1;
            filterOrder = 1;
            % parse removalSteps{~,2} for filtering information
            filterStr = removalSteps{removalStepIdx,2};
            underscores = [strfind(filterStr,'_'),length(filterStr)+1];
            lowIdx = strfind(filterStr,'Low');
            highIdx = strfind(filterStr,'High');
            orderIdx = strfind(filterStr,'Order');
            if(~isempty(lowIdx))
                diff = underscores - lowIdx;
                diff(diff<0) = 10000;
                [~,u] = min(diff);
                u = underscores(u)-1;
                lowPass = str2num(filterStr(lowIdx(1)+3:u));
            end
            if(~isempty(highIdx))
                diff = underscores - highIdx;
                diff(diff<0) = 10000;
                [~,u] = min(diff);
                u = underscores(u)-1;
                highPass = str2num(filterStr(highIdx(1)+4:u));
            end
            if(~isempty(orderIdx))
                diff = underscores - orderIdx;
                diff(diff<0) = 10000;
                [~,u] = min(diff);
                u = underscores(u)-1;
                filterOrder = str2num(filterStr(orderIdx(1)+5:u));
            end
            outputDataTemp = performFilteringStep(outputDataTemp, inputData, highPass, lowPass, filterOrder, folderpath);
        elseif(stepOrder(step) == '1') % curve fitting step
            outputDataTemp = performCurveFittingStep(outputData);
        elseif(stepOrder(step) == '2') % pca step
            outputDataTemp = performPCAStep(outputData);
        end
    end
    
    % store outputDataTemp in outputDataFiltered
    outputDataFiltered{removalStepIdx} = outputDataTemp;
end
warning('on')



end