function [ outputDataFiltered ] = filterArtifactData( outputData, varargin )
% applies a filter to the artifact data in outputData.artifactData, then
% returns the filtered outputData structure

numZeros = 200;

for i = 1:2:length(varargin)
    switch varargin{i}
        case 'filter'
            filterStruct = varargin{i+1};
            b = filterStruct.b;
            a = filterStruct.a;
    end
end
artifactData = outputData.artifactData;
outputDataFiltered = outputData;

% apply filter to all electrodes
for art = 1:numel(artifactData)
    for i = 1:numel(artifactData(art).electrodeNames)
        % reverse the signals
        artifactFlipped = fliplr(squeeze(artifactData(art).artifact(i,:,:)));
        % pad zeros to beginning
        artifactFlipped = [zeros(size(artifactFlipped,1),numZeros)+repmat(mean(artifactFlipped(:,1:10),2),1,numZeros),...
            artifactFlipped];
        % filter and unreverse signal
        filteredSignal = fliplr(filter(b,a,artifactFlipped')');
        outputDataFiltered.artifactData(art).artifact(i,:,:) = filteredSignal(:,1:end-numZeros);
    end
end

end
