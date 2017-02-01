function [ out ] = computeFiringFrequency( cds, stimState,stimDuration, noStimDuration, GTOstim )
% computes the frequency of firing for both the on case (2) and off case
% (1) and outputs this in an nx2 array. Firing frequency is computed as
% number of spikes/time

out = zeros(size(cds.units,2),2);
spikeLabel = zeros(length(cds.units(1).spikes.ts),1);
%% count number of spikes during no stim and stim conditions
if(~GTOstim)
    % find stimulation times for when the stimState changes
    ssDiff = diff(stimState')';
    idx = find(ssDiff);
    timeStim = cds.analog{1,1}.t(ssDiff(:)~=0);
    % count spikes for each time period
    unitCount = 0;
    for unit = 1:size(cds.units,2) % for each unit
        timeStimCount = 1;
        stimOn = 0;
        if(cds.units(unit).ID ~= 255 && cds.units(unit).ID~=0)
            unitCount = unitCount + 1;
            for s = 1:size(cds.units(unit).spikes.ts,1) % for each spike
                % get time of spike
                ts = cds.units(unit).spikes.ts(s);
                % compare ts to next time in timeStim
                if(timeStimCount <= length(timeStim) && ts > timeStim(timeStimCount))
                    timeStimCount = timeStimCount + 1;
                    stimOn = ~stimOn;
                end
                if(stimOn)
                    out(unit,2) = out(unit,2) + 1;
                    spikeLabel(s) = 1;
                else
                    out(unit,1) = out(unit,1) + 1;
                    spikeLabel(s) = 0;
                end
            end
        end
    end
    % stimState has values at the times above
elseif(GTOstim)
% find stimulation times for when the stimState changes
    ssDiff = diff(stimState')';
    idx = find(ssDiff);
    timeStim = cds.lfp.t(ssDiff(:)~=0);
    % count spikes for each time period
    for unit = 1:size(cds.units,2) % for each unit
        if(cds.units(unit).ID ~= 255 && cds.units(unit).ID~=0)
            timeStimCount = 1;
            stimOn = 0;
            for s = 1:size(cds.units(unit).spikes.ts,1) % for each spike
                % get time of spike
                ts = cds.units(unit).spikes.ts(s);
                % compare ts to next time in timeStim
                if(timeStimCount <= length(timeStim) && ts > timeStim(timeStimCount))
                    timeStimCount = timeStimCount + 1;
                    stimOn = ~stimOn;
                end
                if(stimOn)
                    out(unit,2) = out(unit,2) + 1;
                    spikeLabel(s) = 1;
                else
                    out(unit,1) = out(unit,1) + 1;
                    spikeLabel(s) = 0;
                end
            end
        end
    end
end
% figure()
% hold on
% plot(cds.units(1).spikes.ts,spikeLabel,'o')
% plot(cds.analog{1,1}.t,stimState)
out = [out(:,1)./noStimDuration, out(:,2)./stimDuration];

end

