%%%% Raster plot of all spikes for the trial

sortThreshold = 1.0;

%% get data into cds format
funcFolder = pwd;
cds = commonDataStructure();
filepath = 'D:\Lab\Data\SensorStim\Chips_20151123\';
% filename = 'Chips_20151123_TricepsSweep_009'; % NN: 36
% filename = 'Chips_20151123_WristExtVibe_007'; % NN: 3, 93, 4
% filename = 'Chips_20151123_WristFlexVibe_006'; % NN: 21, 35

filename = 'Chips_20151123_GTOStim_03mA_artefactrejection_002';

% filepath = 'D:\Lab\proc-joe-master\SensorStim\Data\Han_20170106\';
% filename = 'Han_20170106_SpindleStim_FCR_area2EMG_003'; % 156
% filename = 'Han_20170106_SpindleStim_FCRdistal_area2EMG_004'; % 59, 109

GTOstim = ~isempty(strfind(lower(filename),'gtostim'));

timeAfterGTOStim = 0.5; % in seconds

labnum = 6;
monkey = 'monkeyHan';
ranBy = 'ranByRaeed';
array = 'arrayLeftS1Area2';
task = 'taskRW';
cds.file2cds([filepath filename], labnum, monkey, task, ranBy, array)

cd(funcFolder);

%% detect when stim is on/off
[stimState,stimDuration,noStimDuration] = determineStimTiming(cds, GTOstim, timeAfterGTOStim);

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
[vals, order] = sort(firingFrequencyRatio(:,1),'descend');
sortFFR = [vals, order];
% clean this up by removing NaN values
sortFFR = sortFFR(~isnan(sortFFR(:,1)),:);

ms=4;
rpLineLength = 0.25;
verticalLine = 1;

figure();
ax1=subplot(2,1,1);
hold on
counter = 1;
for i = sortFFR(:,2)'
    if(sortFFR(counter,1) > sortThreshold)
        if(verticalLine)
            plot([cds.units(i).spikes.ts,cds.units(i).spikes.ts]', ...
            [(counter)*ones(length(cds.units(i).spikes.ts),1)-rpLineLength/2,...
            (counter)*ones(length(cds.units(i).spikes.ts),1)+rpLineLength/2]','k')
        else
            plot(cds.units(i).spikes.ts,(counter)*ones(length(cds.units(i).spikes.ts),1),'k.','markersize',ms)
        end
    end
    counter = counter+1;
end
ylabel('Neuron Number')
% if(GTOstim)
%     title('GTO Stim')
% else
%     title('Spindle Stim')    
% end
title('FCR Han')
% plots spindle or GTO stim state as well
ax2=subplot(2,1,2);
hold on
if(~GTOstim)
    table = cds.analog{1,1};
    plot(table{:,1},table{:,2}/100,'r','linewidth',2)
end
if(GTOstim)
   plot(cds.lfp.t,stimState,'r')
end
xlabel('Time (sec)')
linkaxes([ax1,ax2],'x')