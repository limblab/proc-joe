function [ outputDataFiltered ] = filterArtifactData( outputData, varargin )
% applies a filter to the artifact data in outputData.artifactData, then
% returns the filtered outputData structure

userFilter = 0;
numZeros = 75;

for i = 1:2:length(varargin)
    switch varargin{i}
        case 'filter'
            filterStruct = varargin{i+1};
            userFilter = filterStruct.userFilter;
            if(userFilter)
                b = filterStruct.b;
                a = filterStruct.a;
            end
    end
end
artifactData = outputData.artifactData;
outputDataFiltered = outputData;

if(~userFilter) % use the one here
    sampRate = 30000; % hz
    fc = 1000; % hz
    [b,a] = butter(2,fc/(sampRate/2),'low');
end

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
