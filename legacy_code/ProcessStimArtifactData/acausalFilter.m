function [f] = acausalFilter(data)
    % data is an MxN matrix, M = data, N = number of artifacts (channels,
    % repeated stimuli, etc.). Data has to be along the first dimension
    
    % make filter
    [b,a] = butter(6,[500]/(30000/2),'high');
    % pad data
    data = [zeros(200,size(data,2));data;zeros(200,size(data,2))];
    % acausal filter
    f = fliplr(filter(b,a,fliplr(data')')')';
    
%     [b,a] = butter(2,[7500]/(30000/2),'low');
%     f = fliplr(filter(b,a,fliplr(f')')')';

    % remove padding
    f = f(201:end-200,:);
end

