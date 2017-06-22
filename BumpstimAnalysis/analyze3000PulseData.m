% %% create cds and save cds
% folderpath = 'D:\Lab\Data\StimArtifact\Han\20170621_3000pulses\';
% pwd = cd;
% cd(folderpath);
% fileList = dir('*-s.NEV*');
% filename = fileList(1).name;
% labnum = 6;
% monkey = 'monkeyHan';
% ranBy = 'ranByJoseph';
% array = 'arrayLeftS1';
% task = 'taskRW';
% cds = commonDataStructure();
% cds.file2cds([folderpath,filename],ranBy,array,monkey,labnum,'ignoreJumps',task);
% cd(folderpath)
% save(strcat(filename(1:end-6),'_cds'),'cds','-v7.3');
% cd(pwd);
% 
% load cds if possible
folderpath = 'D:\Lab\Data\StimArtifact\Han\20170621_3000pulses\';
pwd=cd;
cd(folderpath)
fileList = dir('*positiveThreshold_cds.mat');
load(fileList(1).name);
cd(pwd);
% 
%% % break spikes up into each file based on stim time
cd(folderpath)
folderList = dir('chan*');
cd(pwd)
for f = 1:numel(folderList)
%     we used 1000*(f-1) in the times to encode
    data(f).units = cds.units;
    data(f).filename=folderList(f).name;
    for nn = 1:numel(data(f).units)
        spikeMask = data(f).units(nn).spikes.ts > 1000*(f-1) & data(f).units(nn).spikes.ts < 1000*(f);
        data(f).units(nn).spikes = data(f).units(nn).spikes(spikeMask,:);
        data(f).units(nn).spikes.ts = data(f).units(nn).spikes.ts - 1000*(f-1);
    end
end
% 
% %% we need to load in the artifact data
% load(strcat(folderpath,folderList(1).name,filesep,'artifactData.mat'));
% load(strcat(folderpath,folderList(1).name,filesep,'chList.mat'));
% load(strcat(folderpath,folderList(1).name,filesep,'eList.mat'));
% load(strcat(folderpath,folderList(1).name,filesep,'posList.mat'));
% outputData.artifactData = artifactData; clear artifactData;
% outputData.chList = chList; clear chList;
% outputData.eList = eList; clear eList;
% outputData.posList = posList; clear posList;
% 
% % things for Lee
% all artifact shapes -- here
% 
% for fold = 1:numel(folderList)
%     for XLimVar = 1:2
%         disp(fold)
%         load(strcat(folderpath,folderList(fold).name,filesep,'artifactData.mat'));
%         load(strcat(folderpath,folderList(fold).name,filesep,'chList.mat'));
%         load(strcat(folderpath,folderList(fold).name,filesep,'eList.mat'));
%         load(strcat(folderpath,folderList(fold).name,filesep,'posList.mat'));
%         outputData.artifactData = artifactData; clear artifactData;
%         outputData.chList = chList; clear chList;
%         outputData.eList = eList; clear eList;
%         outputData.posList = posList; clear posList;
%         [b,a] = butter(6,500/(30000/2),'high');
%         artifactsPerPlot = 5;
%         numPlots = 1;
%         artifactXData = ((1:1:(size(outputData.artifactData.artifact,3))))/30 - 5;
%         if(XLimVar == 1)
%             figDir = strcat('D:\Lab\Data\StimArtifact\Han\20170621_3000pulses\figures\',strcat(folderList(fold).name,'_short'),filesep);
%             mkdir(figDir(1:end-1));
%             XLim = [-0.1,10];
%         else
%             figDir = strcat('D:\Lab\Data\StimArtifact\Han\20170621_3000pulses\figures\',strcat(folderList(fold).name,'_long'),filesep);
%             mkdir(figDir(1:end-1));
%             XLim = [-0.1,90];
%         end
%         figCounter = 1;
%         for ch = 1:96
%             artIdx = datasample(1:size(outputData.artifactData.artifact,1),artifactsPerPlot*numPlots,'Replace',false);
%             figure('Visible','off')
%             figPrefix = strcat('chan',num2str(ch),'_');
%             for art=1:artifactsPerPlot*numPlots
%                 subplot(2,2,1)
%                 plot(artifactXData,squeeze(outputData.artifactData.artifact(ch,artIdx(art),:)));
%                 hold on
%                 ylabel('Voltage (\muV)')
%                 xlabel('Time (ms)')
%                 title('Raw')
%                 xlim(XLim)
% 
%                 subplot(2,2,2)
%                 plot(artifactXData,squeeze(outputData.artifactData.artifact(ch,artIdx(art),:)));
%                 hold on
%                 ylabel('Voltage (\muV)')
%                 xlabel('Time (ms)')
%                 title('Raw')
%                 xlim(XLim)
%                 ylim([-300,300])
% 
%                 subplot(2,2,3)
%                 plot(artifactXData,fliplr(filter(b,a,fliplr(squeeze(outputData.artifactData.artifact(ch,artIdx(art),:))')))');
%                 hold on
%                 ylabel('Voltage (\muV)')
%                 xlabel('Time (ms)')
%                 title('Filtered')
%                 xlim(XLim)
%                 ylim([-500,500])
% 
%                 subplot(2,2,4)
%                 plot(artifactXData,fliplr(filter(b,a,fliplr(squeeze(outputData.artifactData.artifact(ch,artIdx(art),:))')))');
%                 hold on
%                 ylabel('Voltage (\muV)')
%                 xlabel('Time (ms)')
%                 title('Filtered')
%                 xlim(XLim)
%                 ylim([-300,300])
% 
%                 if(mod(art,artifactsPerPlot)==0)
%                     figName = strcat(figPrefix,num2str(figCounter));
%                     saveFigure(gcf,figDir,figName);
%                     figure('Visible','off')
%                     figCounter = figCounter+1;
%                 end
%             end
%         end
%         close all
%     end
% end

%% raster with stim on as zero for each condition -- here
% fold = 3;
% load(strcat(folderpath,folderList(fold).name,filesep,'artifactData.mat'));
% load(strcat(folderpath,folderList(fold).name,filesep,'chList.mat'));
% load(strcat(folderpath,folderList(fold).name,filesep,'eList.mat'));
% load(strcat(folderpath,folderList(fold).name,filesep,'posList.mat'));
% outputData.artifactData = artifactData; clear artifactData;
% outputData.chList = chList; clear chList;
% outputData.eList = eList; clear eList;
% outputData.posList = posList; clear posList;

%%
nn = 106;
preOffset = 40;
postOffset = 60;    
    
spikes = data(fold).units(nn).spikes.ts;
spikeTimes = zeros(numel(spikes),1);
spikeIdx = zeros(numel(spikes),1);
spikeWaves = zeros(numel(spikes),preOffset+postOffset+1);
mask = zeros(numel(spikes),1);
chan = data(fold).units(nn).chan;
for s = 1:numel(spikes)
    % find stim on time that is appropriate
    stIdx = find(spikes(s) >= outputData.artifactData.stimOn/30000);
    stIdx = stIdx(end);
    spikeIdx(s) = stIdx;
    spikeTimes(s) = spikes(s) - outputData.artifactData.stimOn(stIdx)/30000;
    waveIdx = floor(spikeTimes(s)*30000);
    if(stIdx > 0 && stIdx <= size(outputData.artifactData.stimOn,1) && waveIdx-preOffset >=2 && waveIdx+postOffset <= size(outputData.artifactData.artifact,3)-4)
        spikeWaves(s,:) = squeeze(outputData.artifactData.artifact(chan,stIdx,waveIdx-preOffset:waveIdx+postOffset));
        mask(s,1) = 1;
    else
        spikeWaves(s,:) = -1;
        mask(s,1) = 0;
    end
end    


spikeTimes = spikeTimes*1000;
% 
% plot(spikeTimes(mask==1),spikeIdx(mask==1),'.')
bE = 0:0.5:90;
[bC,~] = histcounts(spikeTimes(mask==1),bE);
bar(bE(1:end-1),bC)
%% plot all waves (from above)

