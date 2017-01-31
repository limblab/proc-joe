function [ stimState,stimDuration,noStimDuration ] = determineStimTiming( cds, GTOstim, timeAfterGTOStim )
%% this function determines when the stimulus is on and stores that in the boolean array stimState

if(~GTOstim)
    timeStim = cds.analog{1,1}.t;
    table = cds.analog{1,1};
    spindleStimValue = table{:,2};
    stimState = zeros(length(timeStim),1);
    numPoints = 100;
    stdStim = std(spindleStimValue(100:100+numPoints));
    pointsLook = 10;
    for i = pointsLook:length(spindleStimValue)
        if(std(spindleStimValue(i-pointsLook+1:i)) > 250*stdStim && max(abs(spindleStimValue(i-pointsLook+1:i))) > 500)
            stimState(i) = 1;
        else
            stimState(i) = 0;
        end
    end
%     figure()
%     hold on
%     plot(spindleStimValue/1000)
%     plot(stimState)
end
if(GTOstim)
    meanNumPoints = 100;

    stimState = zeros(length(cds.lfp.t),1);
    table = cds.lfp;
    meanStim = mean(table{1:meanNumPoints,2});
    stdStim = std(table{1:meanNumPoints,2});
    stimState(:) = abs(meanStim - table{:,2}) > 10*stdStim;
    flagStim = 0;
    countStim = 0;
    timeDiff = cds.lfp.t(2) - cds.lfp.t(1);
    pointsToKeep = timeAfterGTOStim/timeDiff;
    for i = 1:size(cds.lfp.t,1)
        if(stimState(i) == 1) % stim on
            flagStim = 1;
            countStim = 0;
            if(i~=size(cds.lfp.t,1) && stimState(i+1) == 0)
                flagStim = 2;
            end
        elseif(flagStim == 2 && countStim < pointsToKeep) % in zone just after stimulation 
            countStim = countStim + 1;
            stimState(i) = 1;
        else % no stim zone
            flagStim = 0;
            countStim = 0;
        end
    end
    
%     figure();
%     hold on
%     plot(cds.lfp.t,cds.lfp.S1ElecStimSync)
%     plot(cds.lfp.t,stimState*2000)
end

%% compute time of stim and time of no stim
stimDuration = 0;
noStimDuration = 0;

if(~GTOstim)
    timeDiff = cds.analog{1,1}.t(2) - cds.analog{1,1}.t(1);

    for i = 1:length(stimState)
        if(stimState(i) == 1)
            stimDuration = stimDuration + timeDiff;
        else
            noStimDuration = noStimDuration + timeDiff;
        end  
    end
end
if(GTOstim)
    flagStim = stimState(1);
    timeDiff = cds.lfp.t(2) - cds.lfp.t(1);
    for i = 1:size(cds.lfp.t,1)
        if(stimState(i) == 1)
            stimDuration = stimDuration + timeDiff;
        else
            noStimDuration = noStimDuration + timeDiff;
        end
    end 
end
end

