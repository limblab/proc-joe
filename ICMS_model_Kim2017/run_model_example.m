%% define params
PARAMS.gain = 0.2;
PARAMS.spread = 0.25;
PARAMS.refractory_time_constant = 40; % ms
PARAMS.conductance = 2.32;
PARAMS.chronaxie = 0.43; % ms
PARAMS.rheobase = 3.71; % uA
PARAMS.exp_time_constant = 112; % ms

CONSTANTS.t_abs = 1; % ms
CONSTANTS.distances = 1:0.1:3; % mm

%% plot prob_spike vs r for the first pulse for different amplitudes
currents = [10:10:100]; % uA
STIM_PARAMS.pulse_width = 0.2; % ms
STIM_PARAMS.num_pulses = 10; 
STIM_PARAMS.frequency = 200; % Hz

figure()
hold on
colors = [(0+(1:numel(currents))*1/numel(currents))',zeros(numel(currents),1),zeros(numel(currents),1)];
for c = 1:numel(currents)
    STIM_PARAMS.current = currents(c);
    [output_data] = computeNumberSpikes(STIM_PARAMS,PARAMS,CONSTANTS);
    plot(CONSTANTS.distances,output_data.prob_spike(:,1),'color',colors(c,:),'linewidth',1.5)
end

% make legend
legendStr = strsplit(num2str(currents),' ');
for le = 1:numel(legendStr)
    legendStr{le} = [legendStr{le},'\muA'];
end
l=legend(legendStr);
set(l,'box','off','FontSize',14);
xlabel('Distance from stimulated electrode (mm)');
ylabel('Probability of response');
formatForLee(gcf)
set(gca,'fontSize',14)

%% mean spikes as function of distance for different amplitudes
currents = [10:20:90]; % uA
STIM_PARAMS.pulse_width = 0.2; % ms
STIM_PARAMS.num_pulses = 200; 
STIM_PARAMS.frequency = 200; % Hz
num_tests = 100;
total_spikes = zeros(num_tests,numel(currents),numel(CONSTANTS.distances));



for t = 1:num_tests
    for c = 1:numel(currents)
        STIM_PARAMS.current = currents(c);
        [output_data] = computeNumberSpikes(STIM_PARAMS,PARAMS,CONSTANTS);
        total_spikes(t,c,:) = sum(output_data.spike_history,2);
    end
end

mean_total_spikes = squeeze(mean(total_spikes));
std_total_spikes = squeeze(std(total_spikes));
figure()
hold on
colors = [(0+(1:numel(currents))*1/numel(currents))',zeros(numel(currents),1),zeros(numel(currents),1)];
for c = 1:numel(currents)
    plot(CONSTANTS.distances,mean_total_spikes(c,:),'-','color',colors(c,:),'linewidth',1.5)
end
 
% make legend
legendStr = strsplit(num2str(currents),' ');
for le = 1:numel(legendStr)
    legendStr{le} = [legendStr{le},'\muA'];
end
l=legend(legendStr);
set(l,'box','off','FontSize',14);
xlabel('Distance from stimulated electrode (mm)');
ylabel('Mean number of spikes');
formatForLee(gcf)
set(gca,'fontSize',14)

% add error bars
for c = 1:numel(currents)
    for d = 1:numel(CONSTANTS.distances)
        plot([CONSTANTS.distances(d),CONSTANTS.distances(d)],mean_total_spikes(c,d) + [-std_total_spikes(c,d),std_total_spikes(c,d)],'color',colors(c,:),'linewidth',1.5)
    end
end

%% plot effect of varying pulse width on mean spike count for different ICMS amplitudes
currents = [10:10:80]; % uA
pulse_widths = [0.05,0.1,0.2,0.4];

STIM_PARAMS.num_pulses = 200; 
STIM_PARAMS.frequency = 200; % Hz
num_tests = 200;
total_spikes = zeros(num_tests,numel(pulse_widths),numel(currents));



for t = 1:num_tests
    for c = 1:numel(currents)
        for pw = 1:numel(pulse_widths)
            STIM_PARAMS.pulse_width = pulse_widths(pw);
            STIM_PARAMS.current = currents(c);
            [output_data] = computeNumberSpikes(STIM_PARAMS,PARAMS,CONSTANTS);
            total_spikes(t,pw,c) = sum(output_data.spike_history(3,:));
        end
    end
end


mean_total_spikes = squeeze(mean(total_spikes));
std_total_spikes = squeeze(std(total_spikes));
figure()
hold on
colors = [zeros(numel(pulse_widths),1),zeros(numel(pulse_widths),1),(0+(1:numel(pulse_widths))*1/numel(pulse_widths))'];
for pw = 1:numel(pulse_widths)
    plot(currents,mean_total_spikes(pw,:),'-','color',colors(pw,:),'linewidth',1.5)
end

% make legend
legendStr = strsplit(num2str(pulse_widths),' ');
for le = 1:numel(legendStr)
    legendStr{le} = [legendStr{le},'ms'];
end
l=legend(legendStr);
set(l,'box','off','FontSize',14,'location','best');
xlabel('Stimulation amplitude (\muA)');
ylabel('Mean number of spikes');
formatForLee(gcf)
set(gca,'fontSize',14)

% add error bars
for c = 1:numel(currents)
    for pw = 1:numel(pulse_widths)
        plot([currents(c),currents(c)],mean_total_spikes(pw,c) + [-std_total_spikes(pw,c),std_total_spikes(pw,c)],'color',colors(pw,:),'linewidth',1.5)
    end
end

%% plot effect of varying frequency on mean spike count for different ICMS amplitudes
currents = [10:10:80]; % uA
frequencies = [50,100,250,500,1000];

STIM_PARAMS.pulse_width = 0.2;
num_tests = 200;
total_spikes = zeros(num_tests,numel(frequencies),numel(currents));



for t = 1:num_tests
    for c = 1:numel(currents)
        for f = 1:numel(frequencies)
            STIM_PARAMS.frequency = frequencies(f);
            STIM_PARAMS.num_pulses = STIM_PARAMS.frequency;
            STIM_PARAMS.current = currents(c);
            
            [output_data] = computeNumberSpikes(STIM_PARAMS,PARAMS,CONSTANTS);
            total_spikes(t,f,c) = sum(output_data.spike_history(3,:));
        end
    end
end


mean_total_spikes = squeeze(mean(total_spikes));
std_total_spikes = squeeze(std(total_spikes));
figure()
hold on
colors = [zeros(numel(frequencies),1),(0+(1:numel(frequencies))*1/numel(frequencies))',zeros(numel(frequencies),1)];
for f = 1:numel(frequencies)
    plot(currents,mean_total_spikes(f,:),'-','color',colors(f,:),'linewidth',1.5)
end

% make legend
legendStr = strsplit(num2str(frequencies),' ');
for le = 1:numel(legendStr)
    legendStr{le} = [legendStr{le},' Hz'];
end
l=legend(legendStr);
set(l,'box','off','FontSize',14,'location','best');
xlabel('Stimulation amplitude (\muA)');
ylabel('Mean number of spikes');
formatForLee(gcf)
set(gca,'fontSize',14)

% add error bars
for c = 1:numel(currents)
    for f = 1:numel(frequencies)
        plot([currents(c),currents(c)],mean_total_spikes(f,c) + [-std_total_spikes(f,c),std_total_spikes(f,c)],'color',colors(f,:),'linewidth',1.5)
    end
end
