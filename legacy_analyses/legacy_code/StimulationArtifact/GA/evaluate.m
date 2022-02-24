function [ fit ] = evaluate( population )


global N maxBit numPoints minRange maxRange data popSize bFilter aFilter

%% build waveform
waves = zeros(popSize,numPoints);
for i = 1:numPoints
    idxStart = 1+(i-1)*N;
    idxEnd = N*i;
    waves(:,i) = bin2dec(population(:,idxStart:idxEnd))/maxBit*(maxRange-minRange) + minRange;
end

%% filter waveforms
waves = waves - mean(waves,2);
waves = fliplr(filter(bFilter,aFilter,fliplr(waves)')');
%% fitness is sum((yi - y)^2) (distance to each data point, squared)
dataTest = repmat(data,popSize,1);

fit = sum((waves-dataTest).^2,2)/size(data,1);

end

