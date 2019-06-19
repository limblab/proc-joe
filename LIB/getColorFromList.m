function [ color_out ] = getColorFromList( list_idx,color_idx )
% just a function to get a color from a predefined list
% color_idx = 0 is the first entry due to mod weirdness

    switch list_idx
        case 1
            colors_all = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;166,86,40;247,129,191;153,153,153;0,0,0]/255;
            color_out = colors_all(mod(color_idx,size(colors_all,1))+1,:);
        case 2
            colors_all = [230,25,75; 60,180,75; 0,130,200; 245,130,48; 145,30,180; 0 0 0;...
                0 0 128; 240,50,230; 210,245,60; 170,110,40; 0 128,128]/255;
            color_out = colors_all(mod(color_idx,size(colors_all,1))+1,:);
    end
    
    if(color_idx >= size(colors_all,1))
        warning('color idx requested is too large, wrapping around instead')
    end

end

