%% this script looks at the artificial monkey test data. Specifically, it 
% finds spikes at various times after the onset of stimulation. To do this,
% it uses a manual detection of spikes and bursts of spikes, then compares
% spike times (projected based on the artifical monkey's known pattern) to
% the stim times. The data loaded needs to have labels for spike times
% which can be done during initial data processing.

% 1. load in a cds and spike times
folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\artificialMonkey_20171130_dukeBoardV2mod\';
fileNum = 16;

pwd = cd;
cd(folderpath)
fileList = dirSorted('*cds*');

load(fileList(fileNum).name)

%% 
burstIdx = 0:30000*10:30000*60*6;
spikeIdx = 0:30*10:(30*(1000-20));
burstIdx = repmat(burstIdx,numel(spikeIdx),numel(burstIdx));
burstIdx = reshape(burstIdx,size(burstIdx,2)*numel(spikeIdx),1);

spikeIdx = repmat(spikeIdx',numel(burstIdx)/numel(spikeIdx),1);

burstIdxAll = spikeIdx + burstIdx;
%% 2. propagate labels backward and forward in time to get spike times
cds.spikeIdx = cds.spikeBurstTime + burstIdxAll + floor(cds.timeOffset*30000);
cds.spikeIdx = cds.spikeIdx(cds.spikeIdx > 0);

%% save label
save(fileList(fileNum).name,'cds')

%% now we have the spike times for the bursts -- compare to stim onset and 
% find the stim idxs with spikes right after
cds.stimOnIdx = floor(cds.stimOn*30000);
window = [4,5]*30; % in points, not ms, after stim onset

stimIdx = [];

for i = 1:numel(cds.stimOnIdx)
    spikeIdx = cds.spikeIdx - cds.stimOnIdx(i);
    spikeIdx(spikeIdx < min(window)) = 0;
    spikeIdx(spikeIdx > max(window)) = 0;
    if(~isempty(find(spikeIdx~=0)))
        idx = find(spikeIdx~=0);
        stimIdx(end+1,:) = [i, spikeIdx(idx(1))/30];
    end
end

%% plot artifacts 1 at a time for each stimIdx found
saveFigures = 0;

[b,a] = butter(6,500/(30000/2),'high');
xData = ((0:1:size(cds.artifactData.artifact,3)-1) - floor(size(cds.artifactData.artifact,3)/2))/30;
for i = 1:size(stimIdx,1)
    figure
    plot(xData,fliplr(filter(b,a,fliplr(squeeze(cds.artifactData.artifact(stimIdx(i,1),27,:))'))))
    ylim([-400,400])
    xlim([0,14])
    xlabel('Time after stimulation onset (ms)')
    ylabel('Voltage (\muV)')
    formatForLee(gcf)
    
    if(saveFigures)
        figName = strcat('ArtificialMonkey_dukeBoardV2Mod_file',num2str(fileNum),'_window',num2str(window(1)/30),'-',num2str(window(2)/30),'_stim',num2str(stimIdx(i)));
        saveFiguresLIB(gcf,cd,figName);
    end
end
    
%% get non-spike response 
cds.stimOnIdx = floor(cds.stimOn*30000);
window = [4,7]*30; % in points, not ms, after stim onset

stimIdx = [];

for i = 1:numel(cds.stimOnIdx)
    spikeIdx = cds.spikeIdx - cds.stimOnIdx(i);
    spikeIdx(spikeIdx < min(window)) = 0;
    spikeIdx(spikeIdx > max(window)) = 0;
    if(~isempty(find(spikeIdx~=0)))
        idx = find(spikeIdx~=0);
        stimIdx(end+1,:) = [i, spikeIdx(idx(1))/30];
    end
end

responses = fliplr(filter(b,a,fliplr(squeeze(cds.artifactData.artifact(stimIdx(:,1),27,:)))')')';
meanResponse = mean(responses,2);

%% get spike-responses
cds.stimOnIdx = floor(cds.stimOn*30000);
window = [2,4]*30; % in points, not ms, after stim onset

stimIdx = [];

for i = 1:numel(cds.stimOnIdx)
    spikeIdx = cds.spikeIdx - cds.stimOnIdx(i);
    spikeIdx(spikeIdx < min(window)) = 0;
    spikeIdx(spikeIdx > max(window)) = 0;
    if(~isempty(find(spikeIdx~=0)))
        idx = find(spikeIdx~=0);
        stimIdx(end+1,:) = [i, spikeIdx(idx(1))/30];
    end
end
responses = fliplr(filter(b,a,fliplr(squeeze(cds.artifactData.artifact(stimIdx(:,1),27,:)))')')';
responses = responses - meanResponse;
plot(xData,responses(:,[1,3,4,5,6,10]))
ylim([-400,400])
xlim([0,4])

%% look at gain across time after stimulation onset
cds.stimOnIdx = floor(cds.stimOn*30000);
window = [6,10]*30; % in points, not ms, after stim onset

stimIdx = [];

for i = 1:numel(cds.stimOnIdx)
    spikeIdx = cds.spikeIdx - cds.stimOnIdx(i);
    spikeIdx(spikeIdx < min(window)) = 0;
    spikeIdx(spikeIdx > max(window)) = 0;
    if(~isempty(find(spikeIdx~=0)))
        idx = find(spikeIdx~=0);
        stimIdx(end+1,:) = [i, spikeIdx(idx(1))/30];
    end
end

responses = fliplr(filter(b,a,fliplr(squeeze(cds.artifactData.artifact(stimIdx(:,1),27,:)))')')';
%%
xData = ((0:1:size(cds.artifactData.artifact,3)-1) - floor(size(cds.artifactData.artifact,3)/2))/30;

sigAmpMetric = var(responses,1,2);

plotWindow = [find(xData == 1),find(xData == 5)-2];
plot(xData(plotWindow(1):plotWindow(2)),sigAmpMetric(plotWindow(1):plotWindow(2))/max(sigAmpMetric(plotWindow(1):plotWindow(2))),'linewidth',2)
ylabel('Relative gain')
xlabel('Time after stimulation onset (ms)')
formatForLee(gcf)


%% 
window = [6,10]*30; % in points, not ms, after stim onset
responses = [];
fileList = dir('*cds*');

for f = 1:numel(fileList)
    load(fileList(f).name);
    
    cds.stimOnIdx = floor(cds.stimOn*30000);

    stimIdx = [];

    for i = 1:numel(cds.stimOnIdx)
        spikeIdx = cds.spikeIdx - cds.stimOnIdx(i);
        spikeIdx(spikeIdx < min(window)) = 0;
        spikeIdx(spikeIdx > max(window)) = 0;
        if(~isempty(find(spikeIdx~=0)))
            idx = find(spikeIdx~=0);
            stimIdx(end+1,:) = [i, spikeIdx(idx(1))/30];
        end
    end

    responses(:,end+1:end+size(stimIdx,1)) = fliplr(filter(b,a,fliplr(squeeze(cds.artifactData.artifact(stimIdx(:,1),27,:)))')')';
end
%%
