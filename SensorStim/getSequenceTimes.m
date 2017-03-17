function [ sequenceTimes, eventTimes ] = getSequenceTimes(cds,stimState, GTOstim, useEndAsZero )
% this function generates the bins and values for all spikes in cds for a pst
% according to stimState. The outputs are a matrix of bin heights and a 
% matrix of bin widths. THIS IS DONE FOR A GIVEN UNIT NUMBER

%% generate sequences by finding when sequences start and end based on stim state
sequenceTimes = []; % (:,1)=starts, (:,2)=ends, will be adjusted later
eventTimes = []; % (:) = time of stimulation starts
for i = 2:length(stimState)
    if(stimState(i)~=stimState(i-1)) % change in stim state
        % if change is on, add time to sequence start, else add to sequence
        % end
        if(stimState(i) == 1) 
            if(GTOstim)
                sequenceTimes = [sequenceTimes; cds.lfp.t(i), -1];
                if(~useEndAsZero)
                    eventTimes = [eventTimes; cds.lfp.t(i)];
                end
            else
                sequenceTimes = [sequenceTimes; cds.analog{1,1}.t(i), -1];
                if(~useEndAsZero)
                    eventTimes = [eventTimes; cds.analog{1,1}.t(i)];
                end
            end
        else % change is off
            if(GTOstim)
                sequenceTimes(end,2) = cds.lfp.t(i);
                if(useEndAsZero)
                    eventTimes = [eventTimes; cds.lfp.t(i)];
                end
            else
                sequenceTimes(end,2) = cds.analog{1,1}.t(i);
                if(useEndAsZero)
                    eventTimes = [eventTimes; cds.analog{1,1}.t(i)];
                end
            end
        end
    end
end
%% clean times (need to be at least 0.5s apart and at most 1.5s apart)
timeThresh = 0.5; % seconds
timeThresh2 = 3;
i = 2;
while(i < length(eventTimes))
    if(eventTimes(i)-eventTimes(i-1) < timeThresh)
        % remove, don't update i.
        if(i==length(eventTimes))
            eventTimes = eventTimes(1:end-1);
            sequenceTimes = sequenceTimes(1:end-1,:);
        else
            eventTimes = [eventTimes(1:i-1); eventTimes(i+1:end)];
            sequenceTimes = [sequenceTimes(1:i-1,:); sequenceTimes(i+1:end,:)];
            i = i-1;
        end
    elseif(eventTimes(i)-eventTimes(i-1) > timeThresh2)
        % remove i-1, don't update i
        if(i-1 == 1)
            eventTimes = eventTimes(2:end);
            sequenceTimes = sequenceTimes(2:end,:);
            i=i-1;
        else
            % eh do nothing?
        end 
    end
    i = i+1;
end
%% Widen sequences
% add and subtract time from the end and beginning of sequences so that the
% sequences are all as large as possible (same amount of time
% added/subtracted for each)
if(size(sequenceTimes,1) > 1)
    timeShift = max(5,(sequenceTimes(2,1)-sequenceTimes(1,2))/2); 
    dt = timeShift/500;
    flagOverlap = 1;
    while(flagOverlap && timeShift>0)
        sequenceTimesAdjusted = [sequenceTimes(:,1)-timeShift, sequenceTimes(:,2)+timeShift];
        diffTimes = sequenceTimesAdjusted(2:end,1)-sequenceTimesAdjusted(1:end-1,2);
        if(min(diffTimes) > 0)
            flagOverlap = 0;
        else
            flagOverlap = 1;
            timeShift = timeShift-dt;
        end
    end
    timeShift = max(timeShift,0);
    sequenceTimes = [sequenceTimes(:,1)-timeShift, sequenceTimes(:,2)+timeShift];
else
    sequenceTimes = [sequenceTimes(:,1) - 0.8, sequenceTimes(:,2) + 0.8];
end

end

