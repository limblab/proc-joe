function [settlingTime] = computeSettlingTime(data,sampRate,windowSize)
% computes settling time
settlingTime = 0;
meanValue = mean(data(end-windowSize/2:end));
stdValue = std(data(end-windowSize/2:end));

numPoints = 8; % # of points to look at simultaneously
numOutOfRange = 6; % # points out of range that ends the computation 

i = numel(data)-numPoints;
while i > 1 % going backwards through data
    d = data(i:i+numPoints);
    if(sum(d < meanValue-2*stdValue | d > meanValue+2*stdValue) >= numOutOfRange)
        settlingTime = i/sampRate;
        i=0;
    end
    i = i-1;
end

end