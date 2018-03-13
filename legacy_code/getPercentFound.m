function [foundStruct] = getPercentFound(spikeData,artifactData, opts)

    opts = configureOptions(opts);
    
    % 
    foundStruct.percentCorrect = zeros(size(artifactData.artifact,2),numel(unique(artifactData.neuronIdx)));

    for ch = 1:size(artifactData.artifact,2)
        for spike = 1:numel(spikeData.Electrode)
            if(spikeData.Electrode(spike) == ch && spikeData.Unit(spike) == 1)
                foundStruct.percentCorrect(ch,find(artifactData.neuronIdx(spikeData.TimeStamp(spike)) == opts.IDX_ADD_NEURONS)) =  foundStruct.percentCorrect(ch,find(artifactData.neuronIdx(spikeData.TimeStamp(spike)) == opts.IDX_ADD_NEURONS)) + 1;
            end
            
        end
    end
    foundStruct.percentCorrect = foundStruct.percentCorrect/opts.NEURONS_PER_IDX;
    
    
end



function [opts] = configureOptions(optsInput)
    opts.NEURONS_PER_IDX = 300;
    opts.IDX_ADD_NEURONS = [15:60];
    
    %% check if in optsSave and optsSaveInput, overwrite if so
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