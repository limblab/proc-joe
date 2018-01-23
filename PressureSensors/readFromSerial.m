%% reads data from serial port and plots in pseudo realtime
minY = 980;
maxY = 1020;
sampRate = 75;

windowLength = 3; % seconds
xData = ((1:1:sampRate*windowLength)-1)/sampRate;
yData = minY + 10 + zeros(1,numel(xData));

% dataFilename = 'C:\Users\Joseph\Desktop\Lab\Data\PressureSensors\2071207_mineraloiltests\LPS33HW_1.5cmTube_33percentMineralOil_incrementalForce.txt';

s = serial('COM5');
fopen(s);

% fileID = fopen(dataFilename,'w');
try
    flushinput(s);
    linesRead = 0;
    h = figure();
    plotGraph = plot(xData,yData,'k','linewidth',3);
    set(h.Children,'XTick',[]);
    set(h.Children,'YTick',[]);
    set(h.Children,'visible','off');
    ylim([minY,maxY])
    while(ishandle(h))
        line = str2double(fgetl(s));
%         fprintf(fileID,strcat(num2str(line),'\n'));
        yData = circshift(yData,-1);
        yData(end) = line;
        if(~ishandle(h))
            continue;
        end
        set(plotGraph,'YData',yData);
        xData = xData + 1/sampRate;
        drawnow
    end
catch e
    fclose(s);
%     fclose(fileID);
    disp(e.message);
end
    
fclose(s);
% fclose(fileID);

