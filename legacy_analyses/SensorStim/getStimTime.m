function [stimTimes,idxKeep] = getStimTime(cds)
stimTimes = [];

timeStim = cds.analog{1,2}.t;
table = cds.analog{1,2};
spindleStimValue = table{:,2};
stimState = zeros(length(timeStim),1);
numPoints = 100;
maxStim = max(abs(spindleStimValue(100:100+numPoints)));
stdStim = std(spindleStimValue(100:100+numPoints));
pointsLook = 10;
for i = pointsLook:length(spindleStimValue)
    if(max(abs(spindleStimValue(i-pointsLook+1:i))) > maxStim + stdStim*200)
        stimState(i) = 1;
    else
        stimState(i) = 0;
    end
end
timeDiff = cds.analog{1,2}.t(2) - cds.analog{1,2}.t(1);

flagStim = stimState(1);
for i = 1:length(stimState)
    if(stimState(i) == 1 && flagStim == 0) % turns on
        flagStim = 1;
        stimTimes(end+1,1) = i*timeDiff;
    elseif(stimState(i) == 0 && flagStim == 1) % turns off
        flagStim = 0;
        stimTimes(end,2) = i*timeDiff;
    end
end

%% clean times (need to be at least 0.5s apart and at most 1.5s apart)
timeThresh = 0.5; % seconds
timeThresh2 = 3;
i = 1;
idxKeep = 1:1:size(stimTimes,1);
while(i < length(stimTimes)-1)
    if(stimTimes(i,2)-stimTimes(i,1) < timeThresh)
        % join with next in time or just remove, don't update i
        if(stimTimes(i+1,1)-stimTimes(i,1) > 1.5)
            % remove
            
        else
            stimTimes(i,2) = stimTimes(i+1,2);
            if(i==length(stimTimes)-1)
                stimTimes(i,2) = stimTimes(1:end-1,:);
                idxKeep = idxKeep(1:end-1);
            else
                stimTimes = [stimTimes(1:i,:); stimTimes(i+2:end,:)];
                i = i-1;
                idxKeep = [idxKeep(1:i),idxKeep(i+2:end)];
            end 
        end
    end
    i = i+1;
end

end