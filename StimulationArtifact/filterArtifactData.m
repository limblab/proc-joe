function [ outputDataFiltered ] = filterArtifactData( outputData, varargin )
% applies a filter to the artifact data in outputData.artifactData, then
% returns the filtered outputData structure

userFilter = 0;

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
for i = 1:numel(artifactData.electrodeNames)
    % reverse the signals
    artifactFlipped = fliplr(squeeze(artifactData.artifact(i,:,:)));
    % filter and unreverse signal
    outputDataFiltered.artifactData.artifact(i,:,:) = fliplr(filter(b,a,artifactFlipped')');
end


end

