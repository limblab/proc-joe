function [arrayData] = getPValIndependence(arrayData,opts)

    %% configure options and set default values
    opts = configureOpts(opts);
    
     %% for each unit and comparison unit
     for arrayDataIdx = 1:numel(arrayData)
        % initialize independence p vals
        arrayData{arrayDataIdx}.independencePVal = zeros(numel(arrayData),size(arrayData{arrayDataIdx}.bC,1),size(arrayData{arrayDataIdx}.bC,2));
            
        for comparisonIdx = 1:numel(arrayData)
            
            % loop over stim conditions and calculate p vals
            for chan = 1:size(arrayData{arrayDataIdx}.bC,1)
                for wave = 1:size(arrayData{arrayDataIdx}.bC,2)
                    if(arrayDataIdx == comparisonIdx)
                        arrayData{arrayDataIdx}.independencePVal(comparisonIdx,chan,wave) = -100;
                    else
                        arrayData{arrayDataIdx}.independencePVal(comparisonIdx,chan,wave) = ...
                            computeChiSquaredPVal(arrayData{arrayDataIdx}.singleProb.normNumStimsResponsive(chan,wave),arrayData{comparisonIdx}.singleProb.normNumStimsResponsive(chan,wave),...
                                arrayData{arrayDataIdx}.pairwiseProb.bothRespond(comparisonIdx,chan,wave), arrayData{arrayDataIdx}.pairwiseProb.thisRespond(comparisonIdx,chan,wave),...
                                arrayData{arrayDataIdx}.pairwiseProb.comparisonRespond(comparisonIdx,chan,wave),  arrayData{arrayDataIdx}.pairwiseProb.noRespond(comparisonIdx,chan,wave));
                    end
                end
            end
        end
     end
         
end

function [pVal] = computeChiSquaredPVal(singleProb1,singleProb2,bothRespond,respond1,respond2,noRespond)

    Q = (bothRespond-singleProb1*singleProb2)^2/(singleProb1*singleProb2) + ... % both respond
            (respond1-singleProb1*(1-singleProb2))^2/(singleProb1*(1-singleProb2)) + ... % 1 respond
            (respond2-(1-singleProb1)*singleProb2)^2/((1-singleProb1)*singleProb2) + ... % 2 respond
            (noRespond-(1-singleProb1)*(1-singleProb2))/((1-singleProb1)*(1-singleProb2)); % no respond
    pVal = chi2cdf(Q,1,'upper'); % 1 dof

end
function [opts] = configureOpts(optsInput)

    opts = [];
    
    %% check if in opts and optsInput, overwrite if so
    try
        inputFieldnames = fieldnames(optsInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(opts,inputFieldnames{fn}))
               opts.(inputFieldnames{fn}) = optsInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
end