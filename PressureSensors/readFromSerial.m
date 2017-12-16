%% reads data from serial port and plots in pseudo realtime

windowLength = 100;
dataFilename = 'C:\Users\Joseph\Desktop\Lab\Data\PressureSensors\2071207_mineraloiltests\LPS33HW_1.5cmTube_33percentMineralOil_incrementalForce.txt';

yData = zeros(windowLength,1);
s = serial('COM5');
fopen(s);

fileID = fopen(dataFilename,'w');
try
    flushinput(s);
    linesRead = 0;
    h = figure();
    while(ishandle(h))
        line = str2double(fgetl(s));
        fprintf(fileID,strcat(num2str(line),'\n'));
        yData = circshift(yData,-1);
        yData(end) = line;
        if(~ishandle(h))
            continue;
        end
        plot(yData,'k','linewidth',3)
        ylim([950,1300])
        drawnow
    end
catch
    fclose(s);
    fclose(fileID);
    disp('error');
end
    
fclose(s);
fclose(fileID);

