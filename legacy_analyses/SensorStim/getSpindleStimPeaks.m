function [spindlePeakTimes,spindleStimTimes] = getSpindleStimPeaks( cds, stimState, peak )
% calculates the times at which spindle stim is at its peak
% organizes the data based on separate spindle stimulations

spindlePeakTimes = [];
spindleStimTimes = [];
data = cds.analog{1,1};
firstDeriv = diff(data{:,2});
dt = cds.analog{1,1}.t(2) - cds.analog{1,1}.t(1);
waitCounter = 0;
for i = 2:length(stimState)-2
    if(stimState(i) == 1)
        % check for first derivative zero crossing based on which peak to
        % find
        if(peak == 1 && firstDeriv(i-1) > 0 && firstDeriv(i) < 0) % positive peaks
            spindlePeakTimes = [spindlePeakTimes; cds.analog{1,1}.t(i)];
        elseif(peak == -1 && firstDeriv(i-1) < 0 && firstDeriv(i) > 0) % negative peaks
            spindlePeakTimes = [spindlePeakTimes; cds.analog{1,1}.t(i)];
        elseif(peak == 0 && data{i-1,2} < 0 && data{i,2} > 0)
            spindlePeakTimes = [spindlePeakTimes; cds.analog{1,1}.t(i)];
        end
    end
    if(stimState(i-1) == 0 && stimState(i) == 1)
        spindleStimTimes = [spindleStimTimes; cds.analog{1,1}.t(i)];
    end
end
    
    
% 
% figure();
% plot(cds.analog{1,1}.t, cds.analog{1,1}.SpindleStim);
% hold on
% plot(spindlePeakTimes, 2000*ones(length(spindlePeakTimes)),'.','markersize',12);
end

