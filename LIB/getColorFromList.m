function [ color_out ] = getColorFromList( list_idx,color_idx )
% just a function to get a color from a predefined list
% color_idx = 0 is the first entry due to mod weirdness

    switch list_idx
        case 1
            colors_all = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]/255;
            color_out = colors_all(mod(color_idx,numel(colors_all))+1,:);

    end
    
    if(color_idx > numel(colors_all))
        warning('color idx requested is too large, wrapping around instead')
    end

end

