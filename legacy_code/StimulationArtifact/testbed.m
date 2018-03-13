%% salpa -- curve fitting in a window to remove artifact
% the basic idea of this is to window the data, fit a cubic polynomial to
% that window, then use the middle point of that window A(n) to do Y(n) =
% V(n) - A(n). where Y is the output of this, V is the input signal, and A
% is the fit
% 
% this script loads an example artifact, and then tries to remove the
% artifact through the salpa algorithm. The goal of this script is to be a
% test bed for the whole salpa idea

load('artifact_lowGain.mat');
load('Neuron_data_canonical');
artifactData = artifactData(1); 
artifact = squeeze(artifactData.artifact(82,1,100:end)); % start of stimulation to end

% place a spike in the artifact
spikeIdx = 80;
artifact(spikeIdx:spikeIdx+length(neuronMeanWave)-1) = artifact(spikeIdx:spikeIdx+length(neuronMeanWave)-1) + neuronMeanWave'/8;

%% we want to try to do curve fitting on a part
% i have manually said that idx 35 is the start of the exponential decay
startIdx = 3;
endIdx = 20;
data = artifact(startIdx:endIdx);
x = ((startIdx:1:endIdx)');
expEqn = 'poly5';
% expEqn = 'poly3';
f = fit(x,data,expEqn);

y = feval(f,x);
% plot things
figure
subplot(3,1,1)
plot(artifact,'k','linewidth',2)
hold on
plot(x,y,'g','linewidth',1)
subplot(3,1,2)
d=data(:)-feval(f,x);
artSub = artifact;
artSub(startIdx:endIdx) = d;
plot(artSub,'k','linewidth',2)
% acausal filter 
fhigh = 500;
flow = 5000;
[b,a] = butter(6,[fhigh,flow]/(30000/2),'bandpass');
artSub(end:end+200,1) = mean(artSub);
artSubFilt = fliplr(filter(b,a,fliplr(artSub')')')';
subplot(3,1,3)
plot(artSubFilt(1:end-200))
%% now let's try SALPA on the exponential decay aspect
startIdx = 35;
endIdx = 180;
data = artifact(startIdx:endIdx);
x = ((startIdx:1:endIdx)');
y = data;
yFit = [];
expEqn = 'poly3';
N = 50;
for currIdx = 1:numel(data)
    % fit data around currIdx
    windowDataIdx = [currIdx-N/2 currIdx+N/2];
    if(windowDataIdx(1) < 1)
        windowDataIdx = windowDataIdx + abs(windowDataIdx(1)) + 1;
    end
    if(windowDataIdx(2) > numel(data))
        windowDataIdx = windowDataIdx - (windowDataIdx(2) - numel(data));
    end
    windowData = data(windowDataIdx(1):windowDataIdx(2));
    xWindow = x(windowDataIdx(1):windowDataIdx(2));
    
    f=fit(xWindow,windowData,expEqn);
    y(currIdx) = y(currIdx) - feval(f,x(currIdx));
    yFit(end+1) = feval(f,x(currIdx));
end

% plot things
figure
subplot(2,1,1)
plot(artifact)
hold on
plot(x,yFit)
subplot(2,1,2)
plot(x,y)

%% now let's do PCA/ICA, use artifactData.artifact
% grab a single stimulation
load('artifact_lowGain.mat');
load('outputDataExample.mat');
load('Neuron_data_lots');

placeSpikes = 1;
channelsWithSpikes = [1:96]';
stimIdx = 5;
artifactOneStim = squeeze(artifactData.artifact(:,stimIdx,:));

% this is where we would want to place spikes
spikeIdxAll = zeros(size(artifactOneStim,1),1)-1;
waveIdxAll = zeros(size(artifactOneStim,1),1)-1;
if(placeSpikes)
    for ch = channelsWithSpikes'
%         spikeIdx = ceil(rand()*(size(artifactAll,2)-size(neuronMeanWave,2)-2));
%         spikeIdx = ceil(rand()*80);
        spikeIdx = 140;
%         waveIdx = ceil(rand()*size(neuronMeanWave,1));
        waveIdx = 1;
        artifactOneStim(ch,spikeIdx:spikeIdx+size(neuronMeanWave,2)-1) = ... 
            artifactOneStim(ch,spikeIdx:spikeIdx+size(neuronMeanWave,2)-1) + neuronMeanWave(waveIdx,:)/8;
        
        spikeIdxAll(ch) = spikeIdx;
        waveIdxAll(ch) = waveIdx;
    end
end
% run PCA
[coeff,score,latent,tsquared,explained,mu] = pca(artifactOneStim,'Centered',false);
figure
plot(coeff(:,1:1))
coeff(:,1:1) = 0;
artifactPCA = score*coeff';

% plot each channel
eList = outputData.eList;
posList = outputData.posList;
chList = outputData.chList;
figure;
for elec = 1:96
    posIdx=find(strcmp(eList,outputData.artifactData(1).electrodeNames{elec}));
    eRow=posList(posIdx,1);
    eCol=posList(posIdx,2);
    h=subplot(10,10,10*(eRow-1)+eCol);
    hold on
    f=plot(1:1:numel(artifactPCA(elec,:)),artifactPCA(elec,:),'linewidth',2,'color','k');
    if(spikeIdxAll(elec) ~= -1)
        plot(spikeIdxAll(elec):spikeIdxAll(elec)+size(neuronMeanWave,2)-1,artifactPCA(elec,spikeIdxAll(elec):spikeIdxAll(elec)+size(neuronMeanWave,2)-1),'linewidth',2,'color','r')
    end
    ylim([-200,200])
end


%% curve fitting on all channels then PCA with spikes!
load('outputDataExample');
load('Neuron_data_lots');

placeSpikes = 1;
channelsWithSpikes = [33:42,46:49,57:67,72:78,83:86,90:96]';
stimIdx = 5;
artifactOneStim = squeeze(artifactData.artifact(:,stimIdx,100:end));

% this is where we would want to place spikes
spikeIdxAll = zeros(size(artifactOneStim,1),1)-1;
waveIdxAll = zeros(size(artifactOneStim,1),1)-1;
if(placeSpikes)
    for ch = channelsWithSpikes'
%         spikeIdx = ceil(rand()*(size(artifactAll,2)-size(neuronMeanWave,2)-2));
        spikeIdx = ceil(rand()*100);
        waveIdx = ceil(rand()*size(neuronMeanWave,1));
        artifactOneStim(ch,spikeIdx:spikeIdx+size(neuronMeanWave,2)-1) = ... 
            artifactOneStim(ch,spikeIdx:spikeIdx+size(neuronMeanWave,2)-1) + neuronMeanWave(waveIdx,:);
        
        spikeIdxAll(ch) = spikeIdx;
        waveIdxAll(ch) = waveIdx;
    end
end

% curve fitting on all channels
startIdx = 13;
endIdx = 30;
for ch = 1:96
    disp(ch)
    data = artifactOneStim(ch,startIdx:endIdx);
    y=data';
    x = ((startIdx:1:endIdx)');
    expEqn = 'poly4';
    N = 50;
    for currIdx = 1:numel(data)
        % fit data around currIdx
        windowDataIdx = [currIdx-N/2 currIdx+N/2];
        if(windowDataIdx(1) < 1)
            windowDataIdx = windowDataIdx + abs(windowDataIdx(1)) + 1;
        end
        if(windowDataIdx(2) > numel(data))
            windowDataIdx = windowDataIdx - (windowDataIdx(2) - numel(data));
        end
        windowData = data(windowDataIdx(1):windowDataIdx(2));
        xWindow = x(windowDataIdx(1):windowDataIdx(2));

        f=fit(xWindow,windowData',expEqn);
        y(currIdx) = y(currIdx) - feval(f,x(currIdx));
    end
    data = y;
    artifactOneStim(ch,startIdx:endIdx) = data';
end
% plot each channel
eList = outputData.eList;
posList = outputData.posList;
chList = outputData.chList;
figure;
for elec = 1:96
    posIdx=find(strcmp(eList,outputData.artifactData(1).electrodeNames{elec}));
    eRow=posList(posIdx,1);
    eCol=posList(posIdx,2);
    h=subplot(10,10,10*(eRow-1)+eCol);
    hold on
    f=plot(1:1:numel(artifactOneStim(elec,:)),artifactOneStim(elec,:),'linewidth',2,'color','k');
    if(spikeIdxAll(elec) ~= -1)
        plot(spikeIdxAll(elec):spikeIdxAll(elec)+size(neuronMeanWave,2)-1,artifactOneStim(elec,spikeIdxAll(elec):spikeIdxAll(elec)+size(neuronMeanWave,2)-1),'linewidth',2,'color','r')
    end
    ylim([-1000,1000])
end
% 
% % run PCA
% [coeff,score,latent,tsquared,explained,mu] = pca(artifactAll,'Centered',false);
% figure
% plot(coeff(:,1:1))
% coeff(:,1:1) = 0;
% artifactPCA = score*coeff';
% 
% % plot each channel
% eList = outputData.eList;
% posList = outputData.posList;
% chList = outputData.chList;
% fwaves = figure;
% fpower = figure;
% for elec = 1:96
%     posIdx=find(strcmp(eList,outputData.artifactData(1).electrodeNames{elec}));
%     eRow=posList(posIdx,1);
%     eCol=posList(posIdx,2);
%     figure(fwaves)
%     h=subplot(10,10,10*(eRow-1)+eCol);
%     hold on
%     f=plot(1:1:numel(artifactPCA(elec,:)),artifactPCA(elec,:),'linewidth',2,'color','k');
%     if(spikeIdxAll(elec) ~= -1)
%         plot(spikeIdxAll(elec):spikeIdxAll(elec)+size(neuronMeanWave,2)-1,artifactPCA(elec,spikeIdxAll(elec):spikeIdxAll(elec)+size(neuronMeanWave,2)-1),'linewidth',2,'color','r') 
%     end
%     ylim([-1000,1000])
%     
%     figure(fpower)
%     h=subplot(10,10,10*(eRow-1)+eCol);
%     sampRate = 30000;
%     pow=periodogram(artifactPCA(elec,:),[],0:1:2500,sampRate);
%     plot(10*log10(pow),'linewidth',2)
%     xlim([0 2500])
% end

% acausal filter

% plot each channel

%% Ensemble Empirical Mode Decomposition with ICA
[modes,its]=ceemdan(artifact,1,5,10000,1);
[source,A,W] = fastica(modes,'numOfIc',size(modes,1),'displayMode','signals','g','pow3','finetune','tanh','approach','defl','epsilon',0.00001);
%%
figure
subplot(2,2,1)
plot(modes(:,:)')
subplot(2,2,2)
origEmd = A*source;
plot(origEmd(:,:)')
subplot(2,2,3)
origSignal = sum(origEmd,1);
plot(origSignal')
% remove a component from source
subplot(2,2,4)
sourceRemoved = source;
sourceRemoved([2,3,4,5],:) = 0;
origEmdRemoved = A*sourceRemoved;
origSignalRemoved = sum(origEmdRemoved,1);
plot(origSignalRemoved')

%% ICA on all channels across one stimulation
placeSpikes = 1;
channelsWithSpikes = [33:42,46:49,57:67,72:78,83:86,90:96]';
stimIdx = 100;
artifactOneChan = squeeze(artifactData.artifact(73,:,:));
artifactOneStim = squeeze(artifactData.artifact(:,73,:));
% this is where we would want to place spikes
spikeIdxAll = zeros(size(artifactOneStim,1),1)-1;
waveIdxAll = zeros(size(artifactOneStim,1),1)-1;
if(placeSpikes)
    for ch = channelsWithSpikes'
%         spikeIdx = ceil(rand()*(size(artifactAll,2)-size(neuronMeanWave,2)-2));
        spikeIdx = ceil(rand()*80)+40;
        waveIdx = 1;
        artifactOneStim(ch,spikeIdx:spikeIdx+size(neuronMeanWave,2)-1) = ... 
            artifactOneStim(ch,spikeIdx:spikeIdx+size(neuronMeanWave,2)-1) + neuronMeanWave(waveIdx,:);
        
        spikeIdxAll(ch) = spikeIdx;
        waveIdxAll(ch) = waveIdx;
    end
end
%%
[source,A,W] = fastica(artifactOneStim);
plot(source')
%%
sourceMask = min(source,[],2)<-12;
source(sourceMask,:) = 0;
raw = A*source;
% figure
% subplot(2,1,1)
% plot(artifactOneStim(1,:));
% subplot(2,1,2)
% plot(source')

%% denoise wtf
% add noise first?
noise = 1;
figure
decROW = mdwtdec('r',artifact,1,'sym4');
[XD,decDEN] = mswden('den',decROW,'sqtwolog','mln');
resids = artifact-XD;
subplot(3,1,1)
plot(artifact)
subplot(3,1,2)
plot(XD)
subplot(3,1,3)
plot(resids)

%% template subtraction ideas -- cathode for anode, artifact similarity
% should match amplitudes (scale signals) and try to match times -- either
% by doing anode cathode trick or by selecting signals intelligently? then
% apply a filter probably. This is to try to get at neurons during the 
% stimulation. inv(A)*(A(x+n)-y) where Ax = y hopefully. This process would
% leave n untouched while removing the noise x (n is the neuron signal)

load('artifact_lowGain.mat');
load('Neuron_data_canonical');
artifactData = artifactData(1); 
artifact = squeeze(artifactData.artifact(82,1:end,100:end)); % start of stimulation to end
% place a spike in the artifact
% artifact(spikeIdx:spikeIdx+length(neuronMeanWave)-1) = artifact(spikeIdx:spikeIdx+length(neuronMeanWave)-1) + neuronMeanWave'/8;
% plot((artifact-repmat(mean(artifact,2),1,size(artifact,2)))')

anode = artifact(2:2:end,:);
cathode = artifact(1:2:end,:);
spikeIdx = 5;
% amplitude normalization and alignment
cathode = cathode' - repmat(mean(cathode(:,end-50:end),2),1,size(cathode,2))';
anode = anode' - repmat(mean(anode(:,end-50:end),2),1,size(anode,2))';


cathodeSub = cathode;
cathodeSub(spikeIdx:spikeIdx+length(neuronMeanWave)-1,:) = cathodeSub(spikeIdx:spikeIdx+length(neuronMeanWave)-1,:) ... 
    + repmat(neuronMeanWave',1,size(cathodeSub,2));
anodeSub = anode;
% anodeSub(spikeIdx:spikeIdx+length(neuronMeanWave)-1,:) = anodeSub(spikeIdx:spikeIdx+length(neuronMeanWave)-1,:) ... 
%     + repmat(neuronMeanWave,1,size(anodeSub,2));

% for i = 1:size(cathodeSub,2)
%     cathodeSub(:,i) = cathodeSub(:,i).*(cathodeSub(:,i)>0)./repmat((max(cathodeSub(:,i))),size(cathodeSub(:,i),1),1)+...
%         cathodeSub(:,i).*(cathodeSub(:,i)<0)./abs(repmat((min(cathodeSub(:,i))),size(cathodeSub(:,i),1),1));
% end
% for i = 1:size(anodeSub,2)
%     anodeSub(:,i) = anodeSub(:,i).*(anodeSub(:,i)>0)./repmat((max(anodeSub(:,i))),size(anodeSub(:,i),1),1)+...
%         anodeSub(:,i).*(anodeSub(:,i)<0)./abs(repmat((min(anodeSub(:,i))),size(anodeSub(:,i),1),1));
% end

lags = [];

% align cathodeSub and anodeSub based on max point
for i = 2:size(cathodeSub,2)
    [r,lag] = xcorr(cathodeSub(:,1),cathodeSub(:,i));
    [~,idxLag] = max(r);
    % shift cathodeSub(:,i) by idxLag
    cathodeSub(:,i) = circshift(cathodeSub(:,i),lag(idxLag));
    lags(i,1) = lag(idxLag);
end

for i = 2:size(anodeSub,2)
    [r,lag] = xcorr(anodeSub(:,1),anodeSub(:,i));
    [~,idxLag] = max(r);
    % shift cathodeSub(:,i) by idxLag
    anodeSub(:,i) = circshift(anodeSub(:,i),lag(idxLag));
    lags(i,2) = lag(idxLag);
end

% figure;
% subplot(2,1,1)
% plot(cathodeSub(1:size(cathodeSub,2),lags(:,1)==0),'k')
% hold on
% plot(cathodeSub(1:size(cathodeSub,2),lags(:,1)==1),'r')
% 
% subplot(2,1,2)
% plot(anodeSub(:,lags(1:size(anodeSub,2),2)==0),'k')
% hold on
% plot(anodeSub(:,lags(1:size(anodeSub,2),2)==1),'r')

% %
figure;
subplot(3,2,1)
plot(cathodeSub)
subplot(3,2,2)
plot(anodeSub)

templateCathode = -1*mean(anodeSub,2);
templateAnode = -1*mean(cathodeSub,2);

cathodeSub(1:20,:) = cathodeSub(1:20,:) - repmat(templateCathode(1:20,1),1,size(cathodeSub,2));
anodeSub(1:20,:) = anodeSub(1:20,:) - repmat(templateAnode(1:20,1),1,size(anodeSub,2));
% anodeSub = anode - repmat(templateAnode,size(anode,1),1);


subplot(3,2,3)
plot(cathodeSub)
subplot(3,2,4)
plot(anodeSub)

[b,a] = butter(6,[250,10000]/(30000/2),'bandpass');
cathodeFilt = fliplr(cathodeSub(:,:)')';
cathodeFilt = filter(b,a,cathodeFilt);
cathodeFilt = fliplr(cathodeFilt')';
anodeFilt = fliplr(anodeSub(:,:)')';
anodeFilt = filter(b,a,anodeFilt);
anodeFilt = fliplr(anodeFilt')';

subplot(3,2,5)
plot(cathodeFilt)
subplot(3,2,6)
plot(anodeFilt)
hold on
plot(templateCathode,'k','linewidth',2)


%% try longer time filtering idea?
% load once because it takes awhile and im lazy
load('artifactRawExample.mat')
load('Neuron_data_canonical');
elec = 73;
elecName = strcat('PMDelec',num2str(elec));
sampleFreq = 30000;
%% data in cds.lfp{:,2:end}; cds.lfp.(elecName)(.), cds.lfp.t(.) for time
stimStart = artifactData.stimOn(1);
stimData = cds.lfp.(elecName)(:)*8;
offset = 100000;
lenNoStim = numel(stimData(stimStart-5000-offset:stimStart - 5000));
lenStim = numel(stimData(stimStart-4000:stimStart-4000+offset));
len = min(lenStim,lenNoStim);
% add neural data
offset=140;
for i = 1:numel(artifactData.stimOn)
    stimData(artifactData.stimOn(i)+offset:artifactData.stimOn(i)+offset+47) = ...
        stimData(artifactData.stimOn(i)+offset:artifactData.stimOn(i)+offset+47) + neuronMeanWave';
end

for i = 1:numel(artifactData.stimOn)*10
    idx = ceil(rand()*(numel(stimData)-1000));
    stimData(idx:idx+47) = ...
        stimData(idx:idx+47) + neuronMeanWave';
end


% get fft from data before and after stimStart in 50 ms bins and average
idxStart = 1;
lenBin = 30000*50/1000;
f = sampleFreq*(0:floor(lenBin)/2)/lenBin;
powerNoStimAve = zeros(numel(f),1);
powerStimAve = zeros(numel(f),1);
while idxStart < len % do for a certain predefined length
    idxEnd = idxStart + lenBin - 1; % 50 ms bin
    fftNoStim = fft(stimData(idxStart:idxEnd));
    fftNoStim = fftNoStim(1:floor(lenBin/2)+1);
    fftStim = fft(stimData(stimStart-4000+idxStart:stimStart-4000+idxEnd-1));
    fftStim = fftStim(1:floor(lenBin/2)+1);
    powerNoStim = (abs(fftNoStim/lenBin));
    powerStim = (abs(fftStim/lenBin));
    powerNoStimAve = powerNoStimAve + powerNoStim/len;
    powerStimAve = powerStimAve + powerStim/len;
    idxStart = idxStart + lenBin; % update by 50 ms 
end

figure;
subplot(3,1,1)
plot(f,powerStimAve)
subplot(3,1,2)
plot(f,powerNoStimAve)
subplot(3,1,3)
powerRatio = powerNoStimAve./powerStimAve;
plot(f,powerRatio)


filterFreq = powerRatio;

% apply filterFreq to fftStim and ifft to see what we get
fftStim = fft(stimData(stimStart-4000:stimStart-4000+len-1));
fStim = sampleFreq*(0:floor(length(fftStim))/2)/length(fftStim);
filterFreqInter = interp1(f,filterFreq,fStim)';

fftStimFilter = filterFreqInter.*fftStim(1:floor(length(fftStim))/2+1);
iStim = abs(ifft(fftStimFilter,len));
figure
ax(1)=subplot(2,1,1);
plot(1/1000*stimData(stimStart-4000:stimStart-4000+len-1))
ax(2)=subplot(2,1,2);
plot(iStim/1000)
linkaxes(ax,'x');
