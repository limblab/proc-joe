%% spectrum analysis of neuron
counter = 1;
filenames = {'Neuron_data_canonical.mat','Neuron_data_longerTail.mat','Neuron_data_weirdBlip.mat'};
figure
for fname = filenames
    load(fname{1}); % load neuron
    numPoints = 200;
    sampRate = 30000;
    xData = (1:1:numPoints)/sampRate;
    yData = zeros(1,numPoints);
    yData(100:100+length(neuronMeanWave)-1) = yData(100:100+length(neuronMeanWave)-1)+neuronMeanWave;
    
    subplot(2,3,counter)
    plot(xData,yData,'linewidth',2)
    xlabel('Time (s)')
    ylabel('Amplitude (\muV)')
    formatForLee(gcf)
    subplot(2,3,counter+3)
    f=periodogram(yData,[],0:1:5000,sampRate);
    plot(10*log10(f),'linewidth',2)
    xlim([0 5000])
    xlabel('Frequency (Hz)')
    ylabel('Power/Frequency (dB/Hz)')
    formatForLee(gcf)
    counter = counter+1;
end

%% do it for the stimulation artefact
figure;
load('artifact_lowGain.mat')
idx = [52, 71, 82];
sampRate = 30000;
counter = 1;
for i = idx
    yData = squeeze(artifactData(1).artifact(i,1,:));
    xData = (1:1:numel(yData))/sampRate;
    subplot(2,3,counter)
    plot(xData*1000,yData,'linewidth',2)
    xlabel('Time (ms)')
    ylabel('Amplitude (\muV)')
    formatForLee(gcf)
    subplot(2,3,counter + 3)
    f=periodogram(yData,[],0:1:5000,sampRate);
    plot(10*log10(f),'linewidth',2)
    xlim([0 5000])
    xlabel('Frequency (Hz)')
    ylabel('Power/Frequency (dB/Hz)')
    formatForLee(gcf)
    counter = counter+1;
end

%% look at how filter affects neuron and the stimulation artefact
counter = 1;
filenames = {'Neuron_data_canonical.mat','Neuron_data_longerTail.mat','Neuron_data_weirdBlip.mat'};
load('artifact_lowGain.mat')
idx = [52, 71, 82];
sampRate = 30000;
figure
for fname = filenames
    load(fname{1}); % load neuron
    numPoints = 200;
    sampRate = 30000;
    xData = (1:1:numPoints)/sampRate;
    yData = zeros(1,numPoints);
    neuronMeanWave = 150*neuronMeanWave/(max(abs(neuronMeanWave)));
    yData(100:100+length(neuronMeanWave)-1) = yData(100:100+length(neuronMeanWave)-1)+neuronMeanWave;
    
    % filter yData
    sampRate = 30000; % hz
    fhigh = 200;
    flow = 3000;
    [b,a] = butter(6,[fhigh,flow]/(sampRate/2),'bandpass');
    yDataFilt1 = fliplr(yData);
    yDataFilt1 = filter(b,a,yDataFilt1);
    yDataFilt1 = fliplr(yDataFilt1);    

    % plot neuron before and after filtering
    subplot(2,3,counter)
    plot(xData,yData,'linewidth',2,'color','k')
    hold on
    plot(xData,yDataFilt1,'linewidth',2,'color','r')

    xlabel('Time (s)')
    ylabel('Amplitude (\muV)')
    formatForLee(gcf)
    
    % do for artefact
    numPad = 200;
    yData = squeeze(artifactData(1).artifact(idx(counter),1,:))';
    xData = ((1:1:numel(yData)))/sampRate;
    yData(end+1:end+numPad) = mean(yData(end-50:end))+rand(1,numPad);
    threshold = -4*rms(yData(1:97)-mean(yData(1:97)));
    yData = yData - mean(yData(1:97));
    % filter yData
    yDataFilt1 = fliplr(yData);
    yDataFilt1 = filter(b,a,yDataFilt1);
    yDataFilt1 = fliplr(yDataFilt1);
    yData = yData(1:numel(xData));
    yDataFilt1 = yDataFilt1(1:numel(xData));
    
    % plot things
    subplot(2,3,counter+3)
    plot(xData,yData,'linewidth',2,'color','k')
    hold on
    plot(xData,yDataFilt1,'linewidth',2,'color','r');
    xlabel('Time (s)')
    ylabel('Amplitude (\muV)')
    xlim([0 6/1000])
    formatForLee(gcf)
    counter = counter+1;
end

%% do template subtraction and look at artefact power spectrum
figure;
load('artifact_highHeadRoom.mat')
load('Neuron_data_canonical.mat');
idx = [73, 71, 82];
sampRate = 30000;
counter = 1;
neuronIdx = 150;
for i = idx(1)
    template = mean(squeeze(artifactData(1).artifact(i,1:2:end,:)));
    yData = squeeze(artifactData(1).artifact(i,5,:)) - template';
%     yData(neuronIdx:neuronIdx+length(neuronMeanWave)-1) = yData(neuronIdx:neuronIdx+length(neuronMeanWave)-1)+neuronMeanWave';
    xData = ((1:1:numel(yData)))/sampRate;
    
    numPad = 200;
    yData(end+1:end+numPad) = mean(yData(end-50:end))+rand(1,numPad);

    fhigh = 200; % hz
    flow = 1300; % hz
    order = 4;
    [b,a] = butter(order,[fhigh,flow]/(sampRate/2),'bandpass');
    yDataFilt1 = fliplr(yData')';
    yDataFilt1 = filter(b,a,yDataFilt1);
    yDataFilt1 = fliplr(yDataFilt1')';
    yDataFilt1 = yDataFilt1(1:numel(xData));

    [b,a] = butter(order,[fhigh]/(sampRate/2),'high');
    yDataFiltHigh = fliplr(yData')';
    yDataFiltHigh = filter(b,a,yDataFiltHigh);
    yDataFiltHigh = fliplr(yDataFiltHigh')'; 
    yDataFiltHigh= yDataFiltHigh(1:numel(xData));

    [b,a] = butter(order,[flow]/(sampRate/2),'low');
    yDataFiltLow = fliplr(yData')';
    yDataFiltLow = filter(b,a,yDataFiltLow);
    yDataFiltLow = fliplr(yDataFiltLow')'; 
    yDataFiltLow = yDataFiltLow(1:numel(xData));
    
    plot(xData,yData(1:numel(xData)),'linewidth',2,'color','k')
    hold on
    plot(xData,yDataFiltHigh,'linewidth',2,'color','b')
    plot(xData,yDataFilt1,'linewidth',2,'color','r')
    xlabel('Time (s)')
    ylabel('Amplitude (\muV)')
    xlim([0 6/1000])
%     ylim([-1000 1000])
    if(counter == 1)
        l=legend('Artifact','Highpass','Bandpass')
        set(l,'box','off')
    end
    formatForLee(gcf)
%     subplot(2,3,counter + 3)
%     upperLim = 15000;
%     f=periodogram(yData,[],0:1:upperLim,sampRate);
%     plot(10*log10(f),'linewidth',2)
%     xlim([0 upperLim])
%     xlabel('Frequency (Hz)')
%     ylabel('Power/Frequency (dB/Hz)')
%     formatForLee(gcf)
    counter = counter+1;
end

%% bandpass filter multiple neurons and get time shift distribution
load('Neuron_data_lots.mat')
numPoints = 200;
sampRate = 30000;
shifts = zeros(size(neuronMeanWave,1),1);
for i = 1:size(neuronMeanWave,1)
    xData = (1:1:numPoints)/sampRate;
    yData = zeros(1,numPoints);
    nmw = neuronMeanWave(i,:);
    nmw = 150*nmw/(max(abs(nmw)));
    yData(100:100+length(nmw)-1) = yData(100:100+length(nmw)-1)+nmw;

    % filter yData
    sampRate = 30000; % hz
    fhigh = 250;
    flow = 1500;
    [b,a] = butter(6,[fhigh,flow]/(sampRate/2),'bandpass');
    yDataFilt1 = fliplr(yData);
    yDataFilt1 = filter(b,a,yDataFilt1);
    yDataFilt1 = fliplr(yDataFilt1);   
    
    % find shift in peak and save
    [m, minUnfiltered] = min(yData);
    [m, minFiltered] = min(yDataFilt1);
    shifts(i) = (minFiltered-minUnfiltered)/sampRate;
end

[bC,bE] = histcounts(shifts);
bE = bE*1000; % we are in ms now
bFig = bar(bE(1:end-1)+(bE(2)-bE(1))/2,bC,'facecolor','r','edgecolor','r');
bFig.Parent.XTick = [bFig.XData(bFig.YData~=0),0];
f = gcf;
f.Position = [862.6000 322.6000 760.0000 420]
xlabel('Shift Time (ms)')
ylabel('Number of Neurons')
formatForLee(gcf);
title('Bandpass 250-1500Hz')
ylim([0 80])
% xlim([-0.14,0])
%% PCA test for removing artefacts
load('artefact_PowerSpectrum.mat')
idx = [73, 71, 82];

for i = idx(1)
    % get data
    data = squeeze(artifactData(1).artifact(i,:,:));
    % compute PCA for data
    coeff = pca(data');
    dataPCA = coeff'*data;
end