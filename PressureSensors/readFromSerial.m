%% reads data from serial port and plots in pseudo realtime

linesToRead = 10000;
windowLength = 100;

yData = zeros(windowLength,1);
s = serial('COM5');
fopen(s);
try
    flushinput(s);
    linesRead = 0;
    h = figure();
    while(ishandle(h))
        line = str2num(fgetl(s));
        yData = circshift(yData,-1);
        yData(end) = line;
        if(~ishandle(h))
            continue;
        end
        plot(yData)
        ylim([950,1300])
        drawnow
    end
catch
    fclose(s);
    disp('error');
end
    
fclose(s);

