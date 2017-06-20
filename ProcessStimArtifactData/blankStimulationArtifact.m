function [ nevData ] = blankStimulationArtifact( cds, nevData )
% blanks stimulation artifact data
if(any(strcmp('preOffset',fieldnames(cds)))==0)
    cds.preOffset = 27;
end
if(any(strcmp('postOffset',fieldnames(cds)))==0)
    cds.postOffset = 20;
end


preTime = cds.preOffset/30000;
stimOffTime = 0.8/1000;
for st = 1:numel(cds.stimOn)
    stimOnTime = cds.stimOn(st);
    
    % check if any spike in nevData occurs within cds.preOffset + some constant
    % after stimOnTime
    
    wavesIdx = find((nevData.ts(:)-stimOnTime > 0 & nevData.ts(:)-preTime < stimOnTime+stimOffTime));
    
    % for those waves, blank from stimOnTime to stimOnTime + 0.8ms
    for w = 1:numel(wavesIdx)
        numPointsToBlank = ((stimOnTime+stimOffTime) - (nevData.ts(wavesIdx(w))-preTime))*30000;
        nevData.waveforms(wavesIdx(w),1:floor(numPointsToBlank)-1) = 0;
    end
end

end

