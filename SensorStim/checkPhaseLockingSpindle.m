function [ phases,pVal,zScore ] = checkPhaseLockingSpindle(cds,nn,eventStartTimes,stimState)
% checks for phase locking with a neuron by computing the phase for all
% spikes during spindle stimulation (based on analog output data) and then
% using a Rayleigh uniformity test

spindleStim = cds.analog{1,2}.SpindleStim;
t = cds.analog{1,2}.t;
dt = t(2)-t(1);
eventStartIdx = floor(eventStartTimes/dt);
eventEndIdx = find(diff(stimState) == -1);
eventEndTimes = eventEndIdx*dt;
phases = [];
for e = 1:numel(eventStartIdx) % for each event
    %% fit sinusoid to spindle stim data
    % find maxs
    ssData = spindleStim(eventStartIdx(e):eventEndIdx(e));
    tData = t(eventStartIdx(e):eventEndIdx(e));
    maxIdxs = [];

    for idx = 3:numel(ssData)-1
        if(ssData(idx-1) < ssData(idx) && ssData(idx) > ssData(idx+1)) % max
            maxIdxs(end+1) = idx;
        end
    end
    sinFreq = 2*pi/(mean(diff(maxIdxs(2:end)))*dt);
    sinOffset = (maxIdxs(2) + eventStartIdx(e))*dt - pi/2;
    
    figure
    plot(t(eventStartIdx(e):eventEndIdx(e)),800*sin(sinFreq*(t(eventStartIdx(e):eventEndIdx(e))-sinOffset)))
    hold on
    plot(t(eventStartIdx(e):eventEndIdx(e)),spindleStim(eventStartIdx(e):eventEndIdx(e)))
    
    % find spikes within eventStartTimes(e) and eventEndTimes(e)
    spikeMask = cds.units(nn).spikes.ts < eventEndTimes(e) & ...
        cds.units(nn).spikes.ts > eventStartTimes(e);
    spikes = cds.units(nn).spikes.ts(spikeMask);
    % find and store phase for all spikes
    phaseTemp = mod(sinFreq*spikes + sinOffset,2*pi);
    phases = [phases;phaseTemp];
end

figure
bE = linspace(0,2*3.14159,40);
[bC,bE] = histcounts(phases,bE);
bar(bE(1:end-1),bC);

% now we have phases -- do rayleigh test for uniformity
[pVal,zScore] = circ_rtest(phases); % note this func is in ExternalFunctions

end

