function [removalStepsName] = buildRemovalAttemptName(removalSteps)
    removalStepsName = removalSteps;
    removalStepsName{end,end+1} = '';
    for removalAttempt = 1:size(removalSteps,1)
        name = '';
        stepOrder = removalSteps{removalAttempt,1};
        for step = 1:length(stepOrder)
            tempName = '';
            if(stepOrder(step) == '0') % filter, stored in idx 2
                tempName = strcat('filter',strrep(removalSteps(removalAttempt,2),'_',''),'_');
            elseif(stepOrder(step) == '1') % curve fit, stored in idx 3
                tempName = strcat('curveFit',removalSteps(removalAttempt,3),'_');
            elseif(stepOrder(step) == '2') % pca, stored in idx 4
                tempName = strcat('pca',removalSteps(removalAttempt,4),'_');
            end
       
            name = strcat(name,tempName);
        end
        if(length(stepOrder) == 0)
            name = {'nothingDone_'};
        end
        removalStepsName{removalAttempt,end} = name{1}(1:end-1); % removes last '_'
    end
end