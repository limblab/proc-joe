function [ outputData ] = addArtificialNeurons( outputData, meanWave, ampWave, waveIdx, numWaves )
% adds artificial waves to the artifact data. The waveform has the shape of meanWave
% and a maximum magnitude of ampWave. waveIdx refers to the index in the data that 
% the wave should be added to. All channels and all stimulations have the
% same wave added to them
[m,minIdx] = min(meanWave);
for a = 1:numel(outputData.artifactData) 
    % set up times for fake waves. This is redundant honestly but might be
    % used later if the wave timing is not always the same
    outputData.artifactData(a).fakeWaveTimes = zeros(numel(outputData.eList),size(outputData.artifactData(a).artifact,numWaves))-1;
    for i = 1:numel(outputData.eList) % for each electrode  
        for j = 1:size(outputData.artifactData(a).artifact,2) % for each stimulation
            for k = 1:numWaves % for each wave to add
                wave = meanWave/max(abs(meanWave))*ampWave;
                outputData.artifactData(a).artifact(i,j,waveIdx(k):waveIdx(k)+length(wave)-1) = squeeze(outputData.artifactData(a).artifact(i,j,waveIdx(k):waveIdx(k)+length(wave)-1))'+wave;
                outputData.artifactData(a).fakeWaveTimes(i,j,1:numel(waveIdx)) = waveIdx + minIdx;
            end
        end
    end
end


end

