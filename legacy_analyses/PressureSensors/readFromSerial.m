%% reads data from serial port and plots in pseudo realtime
minPressure = 900;
maxPressure = 1600;
minForce = 0;
maxForce = 1023;
sampRate = 50;
max_iter = 10;
iter = 1;
newForce = zeros(max_iter,1);
newPressure = zeros(max_iter,1);

writeFile = 1;
filename = 'C:\Users\Joseph\Desktop\Lab\GIT\proc-joe\PressureSensors\Cat_20180208_soleus_lps22HB_musclesensor031_place1again_20ms_10Hz_10V.txt';

windowLength = 3; % seconds
xData = ((1:1:sampRate*windowLength)-1)/sampRate;
yDataPressure = minPressure + 0.5 + zeros(1,numel(xData));
yDataForce = minForce + 0.5 + zeros(1,numel(xData));

s = serial('COM5'); % Device name in apple OS
fopen(s);

if(writeFile)
    fileID = fopen(filename,'w');
end


try
    flushinput(s);
    linesRead = 0;
    h = figure();
    plotPressure = plot(xData,yDataPressure,'k','linewidth',2);
    hold on
    plotForce = plot(xData,yDataForce,'r','linewidth',2);
    set(h.Children,'XTick',[]);
    %set(h.Children,'YTick',[]);
    set(h.Children,'visible','off');
    ylim([0,1])
    while(ishghandle(h))
        if(s.BytesAvailable ~= 0 && ishghandle(h))
            while(s.BytesAvailable > 0 && iter <= max_iter)
                line = fgetl(s);
                if(writeFile)
                    fprintf(fileID,strcat((line),'\n'));
                end
                idx = strfind(line,',');
                newForce(iter) = (str2double(line(1:idx-1))-minForce)/(maxForce-minForce);
                newPressure(iter) = (str2double(line(idx+1:end))-minPressure)/(maxPressure-minPressure);
                iter = iter + 1;
            end
            iter = iter - 1;
            
            yDataPressure = circshift(yDataPressure,-iter);
            yDataForce = circshift(yDataForce,-iter);
            yDataPressure(end-iter+1:end) = newPressure(1:iter);
            yDataForce(end-iter+1:end) = newForce(1:iter);
            
            if(~ishandle(h))
                continue;
            end
            set(plotPressure,'YData',yDataPressure);
            set(plotForce,'YData',yDataForce);
            xData = xData + 1/sampRate;
            drawnow
            iter = 1;
        else
            pause(0.001);
        end
    end
catch e
    fclose(s);
    if(writeFile)
        fclose(fileID);
    end
    disp(e.message);
end
    
fclose(s);
if(writeFile)
    fclose(fileID);
end