function [ color_out ] = getColorFromList( list_idx,color_idx )
% just a function to get a color from a predefined list
% color_idx = 0 is the first entry due to mod weirdness

    switch list_idx
        case 1
            colors_all = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;166,86,40;247,129,191;153,153,153;0,0,0]/255;
            
        case 2
            colors_all = [230,25,75; 60,180,75; 0,130,200; 245,130,48; 145,30,180; 0 0 0;...
                0 0 128; 240,50,230; 210,245,60; 170,110,40; 0 128,128]/255;
            
        case 3 % 3 color red gradient
            colors_all = [253,204,138; 252,141,89; 215,48,31]/255;
        case 4  % 4 color blue gradient
            colors_all = [189,215,231; 107,174,214; 49,130,189; 8,81,156]/255;
        case 5 % 5 colors, inferno gradient
            colors_all = [0.0015    0.0005    0.0139
                        0.2582    0.0386    0.4065
                        0.5783    0.1480    0.4044
                        0.8650    0.3168    0.2261
                        0.9876    0.6453    0.0399
                        0.9884    0.9984    0.6449];
        case 6 % 4 color viridis gradient
            colors_all = [0.2670    0.0049    0.3294
                        0.1906    0.4071    0.5561
                        0.2080    0.7187    0.4729
                        0.9932    0.9062    0.1439];
                    
        case 7 % same as case 1, but lighter
            colors_all = ([228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;166,86,40;247,129,191;153,153,153;0,0,0]/255);
            colors_all = colors_all + (1-colors_all)*1/4; 
            colors_all(colors_all < 0) = 0;
            colors_all(colors_all > 1) = 1;
            
        case 8 % same as case 1, but darker
            colors_all = ([228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;166,86,40;247,129,191;153,153,153;0,0,0]/255);
            colors_all = colors_all*5/8; 
            colors_all(colors_all < 0) = 0;
    end
    color_out = colors_all(mod(color_idx,size(colors_all,1))+1,:);
    if(color_idx >= size(colors_all,1))
        warning('color idx requested is too large, wrapping around instead')
    end

end

