folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20180311_stimswitch\';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';

inputData.array = 'arrayLeftS1';
inputData.monkey = 'monkeyHan';
inputData.ranBy = 'ranByJoe';
inputData.lab = 6;
inputData.mapFile = strcat('mapFile',mapFileName);
inputData.task = 'taskRW';

pwd=cd;
cd(folderpath)
fileList = dir('*nev*');
%%
cds = commonDataStructure();
cds.file2cds(strcat(folderpath,fileList(1).name),inputData.array,inputData.monkey,inputData.ranBy,...
    inputData.lab,inputData.mapFile,inputData.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
cd(pwd);


%%
syncName = 'ainp16';
stimOn=find(diff(cds.analog{1}.(syncName)-mean(cds.analog{1}.(syncName))>3)>.5);
        stimOff=find(diff(cds.analog{1}.(syncName)-mean(cds.analog{1}.(syncName))<-3)>.5);
        
ts = 25;
artifact = zeros(size(cds.lfp,2)-1,numel(stimOn),30*ts); % 5ms of data
for i = 1:numel(stimOn)
%     artifact(1,i,:) = cds.analog{1,1}.ainp15(stimOn(i):stimOn(i)+30*ts-1);
    artifact(:,i,:) = cds.lfp{stimOn(i):stimOn(i)+30*ts-1,2:end}';
end

        
        
%% 
idxPlot = 100;
numPlot = 10;
plotFiltered = 0;
ts = 10;

xData = ((1:size(artifact,3))-1)/30 - 25;
figure();
% subplot(2,1,1)
if(plotFiltered)
    plot(xData,acausalFilter(squeeze(artifact(1,idxPlot:idxPlot+numPlot-1,:))'));
else
%     plot(xData,filtfilt(b,a,squeeze(artifact(1,idxPlot:idxPlot+numPlot-1,:))'));
    plot(xData,squeeze(artifact(1,idxPlot:idxPlot+numPlot-1,:))');
end
ylim([-8000,8000])
xlim([0,min(ts,8)])
formatForLee(gcf)
% subplot(2,1,2)
% 
% 
% if(plotFiltered)
%     plot(xData,acausalFilter(squeeze(artifact(2,idxPlot:idxPlot+numPlot,:))'));
% else
%     plot(xData,squeeze(artifact(2,idxPlot:idxPlot+numPlot,:))');
% end
% ylim([-300,300])
% xlim([0,min(ts,8)])
% 
% formatForLee(gcf)



%% 
t = cds.lfp{:,1};
filteredData = cds.lfp{:,2:end};
syncName = 'ainp16';
stimOn=find(diff(cds.analog{1}.(syncName)-mean(cds.analog{1}.(syncName))>3)>.5);
        stimOff=find(diff(cds.analog{1}.(syncName)-mean(cds.analog{1}.(syncName))<-3)>.5);
        
ts = 10;
        
artifact = zeros(size(filteredData,2),numel(stimOn),30*ts); % 5ms of data
for i = 1:size(filteredData,2)
    filteredData(:,i) = acausalFilter(filteredData(:,i));
    for j = 1:numel(stimOn)
        artifact(i,j,:) = filteredData(stimOn(j):stimOn(j)+30*ts-1,i);
    end
end
xData = (0:1:30*ts-1)/30;
%%
figure
% xData = (0:1:size(artifact,3)-1)/30 - size(artifact,3)/2/30;
xData = (0:1:size(artifact,3)-1)/30;

for i = 1:size(artifact,1)
    subplot(6,6,i)
    plot(xData,squeeze(artifact(i,1:2:10,:)));
    ylim([-10000,10000])
%     ylim([-500,500])
    xlim([0,5])
    title(num2str(i));
    
end
        

%%  extract artifact data if required for whatever reason

% need to manually load waveforms file
artifact = cds.artifactData.artifact;
artifact = artifact(:,11,:);
artifact = permute(artifact,[2,1,3]);

% artifact_wave1 = artifact(:,waveforms.waveSent==1,:);
% artifact_wave2 = artifact(:,waveforms.waveSent==2,:);
% artifact = artifact_wave1;