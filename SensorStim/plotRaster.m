function [ figHandle,sortFFR ] = plotRaster(cds, sortThreshold, GTOstim, timeAfterGTOStim,varargin)
% plots a raster plot given the commondatastructure, a sortThreshold,
% whether this is spindle stim or GTOstim, and time after GTOstim to look
% at
specificNeurons=0;
neurons = [];
for i = 1
    switch varargin{i}
        case 'specificNeurons'
            neurons = varargin{i+1};
            specificNeurons = 1;
    end
end

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


%% plot raster
[vals, order] = sort(firingFrequencyRatio(:,1),'descend');
sortFFR = [vals, order];
% clean this up by removing NaN values
sortFFR = sortFFR(~isnan(sortFFR(:,1)),:);

ms=6;
rpLineLength = 0.25;
verticalLine = 0;

figHandle = figure();
ax1=subplot(2,1,1);
hold on
counter = 1;
for i = sortFFR(:,2)'
    if(specificNeurons)
        if(sum(i==neurons)>0)
            plot(cds.units(i).spikes.ts,(counter)*ones(length(cds.units(i).spikes.ts),1),'k.','markersize',ms)
            counter = counter+1;
        end
        counter = counter-1;
    elseif(sortFFR(counter,1) > sortThreshold)
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
ylabel('Stimulation (V)')
linkaxes([ax1,ax2],'x')


end

