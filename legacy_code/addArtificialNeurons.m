function [ outputData ] = addArtificialNeurons(artifact,neuronData,opts)
% adds artificial waves to the artifact data. The waveform has the shape of meanWave
% and a maximum magnitude of ampWave. waveIdx refers to the index in the data that 
% the wave should be added to. All channels and all stimulations have the
% same wave added to them

    opts = configureOptions(opts);

    [~,minIdx] = min(neuronData.meanWave);
    neuronData.meanWave = neuronData.meanWave/(max(abs(neuronData.meanWave)))*opts.WAVE_AMPLITUDE;
    
    % we need opts.NEURONS_PER_IDX neurons at each idx
    outputData.artifact = zeros(opts.NEURONS_PER_IDX*numel(opts.IDX_ADD_NEURONS),size(artifact,2),size(artifact,3));
    outputData.neuronIdx = zeros(opts.NEURONS_PER_IDX*numel(opts.IDX_ADD_NEURONS),1);

    for i = 1:numel(opts.IDX_ADD_NEURONS)
        idx = opts.IDX_ADD_NEURONS(i) - minIdx + 1;
        % sample artifacts on all channels
        for ch = 1:size(artifact,2)
            data = datasample(artifact(:,ch,:),opts.NEURONS_PER_IDX,'replace',false);
            data(:,1,idx:idx+length(neuronData.meanWave)-1) = data(:,1,idx:idx+length(neuronData.meanWave)-1) + reshape(neuronData.meanWave,1,1,length(neuronData.meanWave));
            outputData.artifact(1+(i-1)*opts.NEURONS_PER_IDX:(i)*opts.NEURONS_PER_IDX,ch,:) = data;
            outputData.neuronIdx(1+(i-1)*opts.NEURONS_PER_IDX:(i)*opts.NEURONS_PER_IDX,1) = opts.IDX_ADD_NEURONS(i);
        end 
    end

end


function [opts] = configureOptions(optsInput)

    opts.IDX_ADD_NEURONS = [20:20:150];
    opts.NEURONS_PER_IDX = 300;
    opts.WAVE_AMPLITUDE = 200;
    
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
