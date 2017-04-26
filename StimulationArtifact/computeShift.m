function [shift] = computeShift(meanWave,filterStr, inputData)

% pad meanWave
padVal = mean(meanWave(1:10));
meanWave = [ones(1,50)*padVal, meanWave];

% split filterName into filter components
lowPass = -1;
highPass = -1;
filterOrder = 1;
% parse filterStr for filtering information
underscores = [strfind(filterStr,'_'),length(filterStr)+1];
lowIdx = strfind(filterStr,'Low');
highIdx = strfind(filterStr,'High');
orderIdx = strfind(filterStr,'Order');
if(~isempty(lowIdx))
    diff = underscores - lowIdx;
    diff(diff<0) = 10000;
    [~,u] = min(diff);
    u = underscores(u)-1;
    lowPass = str2num(filterStr(lowIdx(1)+3:u));
end
if(~isempty(highIdx))
    diff = underscores - highIdx;
    diff(diff<0) = 10000;
    [~,u] = min(diff);
    u = underscores(u)-1;
    highPass = str2num(filterStr(highIdx(1)+4:u));
end
if(~isempty(orderIdx))
    diff = underscores - orderIdx;
    diff(diff<0) = 10000;
    [~,u] = min(diff);
    u = underscores(u)-1;
    filterOrder = str2num(filterStr(orderIdx(1)+5:u));
end

filterType = 'n';
if(highPass == -1 && lowPass == -1)
    filterType = 'n';
elseif(highPass == -1)
    filterType = 'low';
    filterCutoff = lowPass;
elseif(lowPass == -1)
    filterType = 'high';
    filterCutoff = highPass;
else
    filterType = 'bandpass';
    filterCutoff = [highPass,lowPass];
end

if(strcmp(filterType,'n')==1) % check for invalid filter, if so shift = 0 and return
    shift = 0;
else % if not invalid filter, build and compute shift
    
    [b,a] = butter(filterOrder,filterCutoff/(30000/2),filterType);
    filtMeanWave = fliplr(meanWave);
    filtMeanWave = filter(b,a,filtMeanWave);
    filtMeanWave = fliplr(filtMeanWave);
    
    % compute shift in min amplitude idx and put that into shift
    [m,minMeanWaveIdx] = min(meanWave);
    [m,minFiltWaveIdx] = min(filtMeanWave);
    shift = floor(minFiltWaveIdx - minMeanWaveIdx); % treated as shift in index not time everywhere
end


end