function [ neuronsTotal, neuronsFound, thresholdCrossings ] = thresholdData(data,fakeWaveIdx, threshold,  tolerance, shift, inputData )
%
% perform thresholding to get threshold crossings
% threshold = thresholdMult*rms(data(inputData.presample + inputData.windowSize/2:end));
% find all indices with a threshold crossing. Do this by sliding a window backwards,
% filtering with a high pass at like 1 Hz, and then looking at threshold
% crossings based on rms in that window.
neuronsFound = 0;
crossingIdx = [-1];

threshold = rms(data(end-100:end))*-4;
% find threshold based on last (inputData.presample) data points
windowData = data(inputData.presample:end);
for i = 1:numel(windowData)
    if(windowData(i) < threshold)
        crossingIdx(end+1,1) = i+inputData.presample-1;
    end
end

crossingIdx = unique(crossingIdx); % remove duplicates
% remove indexes with long strings of crossings
maxLength = 5;
i = 2;
currentChainLength = 1;
while i < numel(crossingIdx)
    if(crossingIdx(i) == crossingIdx(i-1) + 1)
        currentChainLength = currentChainLength + 1;
    elseif(currentChainLength >= maxLength)
        % remove :D
        endIdx = i;
        startIdx = i-currentChainLength;
        if(startIdx == 1)
            crossingIdx = crossingIdx(endIdx+1:end);
        elseif(endIdx == numel(crossingIdx))
            crossingIdx = crossingIdx(1:startIdx-1);
        else
            crossingIdx = [crossingIdx(1:startIdx-1);crossingIdx(endIdx:end)];
        end
        i = startIdx-1;
        currentChainLength = 1;
    else
        currentChainLength = 1;
    end
    i=i+1;
end
% remove -1 entry
if(numel(crossingIdx) > 1)
    crossingIdx = crossingIdx(2:end);
else
    crossingIdx = [];
end
% remove magnitudes that are too large
i=1;
while i <= numel(crossingIdx)
    if(abs(data(crossingIdx(i))-mean(data(inputData.presample+100:end))) > 100)
        if(i == 1 && numel(crossingIdx) > 1)
            crossingIdx = crossingIdx(2:end);
        elseif(i == 1)
            crossingIdx = [];
        elseif(i == numel(crossingIdx))
            crossingIdx = crossingIdx(1:end-1);
        else
            crossingIdx = [crossingIdx(1:i-1); crossingIdx(i+1:end)];
        end
        i=i-1;
    end
    i = i+1;
end
% check if fakeWaveIdx is in crossingIdx
for i = 1:numel(fakeWaveIdx)
    if(sum(fakeWaveIdx(i) - crossingIdx + shift <= tolerance & fakeWaveIdx(i)-crossingIdx+shift >= -tolerance) > 0)
        neuronsFound = neuronsFound+1;
    end
end

% store some output data
thresholdCrossings = numel(crossingIdx); 
neuronsTotal = numel(fakeWaveIdx);

end

