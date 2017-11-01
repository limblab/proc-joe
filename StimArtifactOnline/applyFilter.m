function [ filteredData ] = applyFilter( data )
% this function takes in the raw data and applies a backwards causal
% filter then outputs the results

%%%%% DATA SHOULD COME IN AS A ROW VECTOR

% create filter
[bFilter,aFilter] = butter(6,500/(30000/2),'high');

% pad data
numPad = 200;
dataTemp = [data(1,:),mean(data(1,end-min(length(data(1,:))-1,20):end))*ones(1,numPad)];

% flip in time, apply filter, flip back in time
dataTemp = fliplr(filter(bFilter,aFilter,fliplr(dataTemp)));

% remove padded values
filteredData(1,:) = dataTemp(1,1:end-numPad);

% filteredData is then returned

end

