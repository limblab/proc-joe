function [f] = acausalFilter(data)
    [b,a] = butter(6,[500]/(30000/2),'high');
    data = [zeros(200,size(data,2));data;zeros(200,size(data,2))];
    f = fliplr(filter(b,a,fliplr(data')')')';
    
    [b,a] = butter(2,[7500]/(30000/2),'low');
    f = fliplr(filter(b,a,fliplr(f')')')';

    f = f(201:end-200,:);
    
end

