function [ density ] = computeDensity(data,cutoff)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
binSize = 1; % its relative, so I chose 1
density = 0;
totalDensity = 0;
for i = 1:numel(data)
    if(data(i) > cutoff)
        density = density + data(i)*binSize;
    end
    totalDensity = totalDensity + data(i)*binSize;
    
end

density = density/totalDensity;

end

