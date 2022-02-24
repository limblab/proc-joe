function [ figHandle,sortFFR ] = plotRaster(cds, neuronNumber, eventTimes, sequenceTimes, varargin)
% plots a raster plot given the commondatastructure, a sortThreshold,
% whether this is spindle stim or GTOstim, and time after GTOstim to look
% at
cRect = [255 140 140]/256;
plotStimTime = 0;
stimTimes = [];
for i = 1:2:size(varargin,2)
    switch varargin{i}
        case 'plotStimTime'
            plotStimTime = 1;
            stimTimes = varargin{i+1};
    end
end

ms=4;
rpLineLength = 0.33;
verticalLine = 1;

figure();
hold on
counter = 1;
preTime = max(eventTimes(:) - sequenceTimes(:,1));
postTime = max(sequenceTimes(:,2)-eventTimes(:));
for i = 1:numel(eventTimes)
    spikeMask = (cds.units(neuronNumber).spikes.ts  > eventTimes(i) - preTime & ...
            cds.units(neuronNumber).spikes.ts < eventTimes(i) + postTime);
    spikesPlot = cds.units(neuronNumber).spikes.ts(spikeMask) - eventTimes(i);
%     plot(spikesPlot*1000,(i)*ones(length(spikesPlot),1),'k.','markersize',ms)
    if(plotStimTime)
        timePlot = (stimTimes(i,2) - stimTimes(i,1))*1000; % s to ms
        hold on
        % plot box from 0 to timePlot, and with height = 1 (0.5 around
        % center)
        pRect = [0, i-0.5, timePlot, 1]; % [x y dx dy]
        if(i==1)
            pRect = [0,i-1,timePlot,1.5];
        end
        r=rectangle('Position',pRect,'FaceColor',cRect,'EdgeColor',cRect);       
    end

    plot([spikesPlot*1000,spikesPlot*1000]', ...
        [(i)*ones(length(spikesPlot),1)-rpLineLength/2,...
        (i)*ones(length(spikesPlot),1)+rpLineLength/2]','k','linewidth',1.5)
    
end
ylim([0,numel(eventTimes)]+rpLineLength/2)
xlim([-1000*preTime,1000*postTime])
ylabel('Stimulation Trial')
xlabel('Time after stimulation onset (ms)')
formatForLee(gcf);
ha = gca;
ha.YRuler.MinorTickValues = [1:1:numel(eventTimes)];
end

