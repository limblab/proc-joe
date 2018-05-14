function [ out ] = isEqual( data1, data2 )

    out = 0;
    if(~isempty(data1) && ~isempty(data2))
        if(mean(data1) == 0)
            off=0.001;
        else
            off = mean(data1(~isnan(data1))).*0.001;
        end
        out = data1 > data2-off & data1 < data2+off;
    end
    
    if(size(out,1) > 1)
        out = sum(out) > 0;
    end
end

