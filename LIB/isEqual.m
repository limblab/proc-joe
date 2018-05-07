function [ out ] = isEqual( data1, data2 )

    out = 0;
    if(~isempty(data1) && ~isempty(data2))
        off = mean(data1).*0.001;
        out = data1 > data2-off & data1 < data2+off;
    end
    
end

