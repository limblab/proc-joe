function [] = generateRaster(cds,neuronNumber,eventTimes,sequenceTimes,varargin)

plotStimTime = 0;
GTOstim = 0;
for i = 1:2:size(varargin,2)
    switch varargin{i}
        case 'plotStimTime'
            plotStimTime = varargin{i+1};
        case 'GTOstim'
            GTOstim = varargin{i+1};
    end
    
end

if(plotStimTime)
    if(~GTOstim)
        [stimTimes,idxKeep] = getStimTime(cds);
%         stimTimes(:,1) = eventTimes(idxKeep');
    else
        stimTimes(:,1) = eventTimes;
        stimTimes(:,2) = eventTimes + 20/1000; % 20 ms
    end
    plotRaster(cds,neuronNumber,eventTimes,sequenceTimes,'plotStimTime',stimTimes)
else
    plotRaster(cds,neuronNumber,eventTimes,sequenceTimes)
end
% call plotRaster

end