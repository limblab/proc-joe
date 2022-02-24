%%%% Raster plot of all spikes for the trial

sortThreshold = 1.0;

%% get data into cds format
funcFolder = pwd;
cds = commonDataStructure();
filepath = 'D:\Lab\Data\SensorStim\Han_20170209\';

filename = 'Han_20170209_GTOstim_ECU_08mA_area2_001.mat';

GTOstim = ~isempty(strfind(lower(filename),'gtostim'));

timeAfterGTOStim = 0.01; % in seconds

labnum = 6;
monkey = 'monkeyHan';
ranBy = 'ranByRaeed';
array = 'arrayLeftS1Area2';
task = 'taskRW';
cds.file2cds([filepath filename], labnum, monkey, task, ranBy, array)

cd(funcFolder);

%% detect when stim is on/off
[stimState,stimDuration,noStimDuration] = determineStimTiming(cds, GTOstim, timeAfterGTOStim);
[sequenceTimes, eventTimes] = getSequenceTimes(cds, stimState,GTOstim,0);
%% compute firing frequency for all units against time during and not during stimulation
% firing frequency = number of spikes/time of stim or number of spikes/time
% of not stim

% computeFiringFrequency returns an nx2 array of firing rates for the off
% case (1) and on case (2)

firingFrequency = computeFiringFrequency(cds, stimState, stimDuration, noStimDuration, GTOstim);
firingCounts = firingFrequency.*repmat([noStimDuration, stimDuration],size(firingFrequency,1),1);
% ratio is during stim divided by not during stim
firingFrequencyRatio = [firingFrequency(:,2)./firingFrequency(:,1), (1:1:length(firingFrequency))'];
%% raster: all units against time sorted by firing frequency ratio
% sort firingFrequencyRatio while keeping indexes preserved
% clean this up by removing NaN values
nn = 73;

ms=4;
rpLineLength = 0.25;
verticalLine = 1;

figure();
hold on
counter = 1;
preTime = max(eventTimes(:) - sequenceTimes(:,1));
postTime = max(sequenceTimes(:,2)-eventTimes(:));
for i = 1:numel(eventTimes)

    spikeMask = (cds.units(nn).spikes.ts  > eventTimes(i) - preTime & ...
            cds.units(nn).spikes.ts < eventTimes(i) + postTime);
    spikesPlot = cds.units(nn).spikes.ts(spikeMask) - eventTimes(i);
%     plot(spikesPlot*1000,(i)*ones(length(spikesPlot),1),'k.','markersize',ms)
    plot([spikesPlot,spikesPlot]', ...
        [(i)*ones(length(spikesPlot),1)-rpLineLength/2,...
        (i)*ones(length(spikesPlot),1)+rpLineLength/2]','k')
     
end
ylim([1-rpLineLength/2,numel(eventTimes)]+rpLineLength/2)
ylabel('Stimulation Trial')
xlabel('Time after stimulation (ms)')
% if(GTOstim)
%     title('GTO Stim')
% else
%     title('Spindle Stim')    
% end
% title('FCR Han')
% % plots spindle or GTO stim state as well
% ax2=subplot(2,1,2);
% hold on
% if(~GTOstim)
%     table = cds.analog{1,1};
%     plot(table{:,1},table{:,2}/100,'r','linewidth',2)
% end
% if(GTOstim)
%    plot(cds.lfp.t,stimState,'r')
% end
% xlabel('Time (sec)')
% linkaxes([ax1,ax2],'x')