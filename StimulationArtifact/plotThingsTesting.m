pwd=cd;
folderpath = 'D:\Lab\proc-joe\StimulationArtifact\';
filename = 'Han_20170516.nev';
labnum = 6;
monkey = 'monkeyHan';
ranBy = 'ranByTucker';
array = 'arrayLeftS1';
task = 'taskRW';
cds=commonDataStructure();
cds.file2cds([folderpath,filename],ranBy,array,monkey,labnum,'ignoreJumps',task);
cd(pwd);
%%
folderpath = 'D:\Lab\Data\StimArtifact\Han\noStim\';
filename = 'Han_20170516_freeReaching-handHolding-forTreats_059_cds.mat';
load([folderpath,filename])

%%
[b,a] = butter(6,500/30000*2,'high');
numZeros = 200;
data = cds.units(38).spikes{:,2:end};

artifactFlipped = fliplr(data);
artifactFlipped = [zeros(size(artifactFlipped,1),numZeros)+repmat(mean(artifactFlipped(:,1:10),2),1,numZeros),...
    artifactFlipped,zeros(size(artifactFlipped,1),numZeros)+repmat(mean(artifactFlipped(:,end-10:end),2),1,numZeros)];
dataFilt = fliplr(filter(b,a,artifactFlipped')');
plotX = (0:1:47)/30000*1000;
figure
plot(plotX,data(1:10,:)'*8)
figure
plot(plotX,dataFilt(1:10,numZeros-5:numZeros+42)'*8)
xlabel('Time (ms)')
ylabel('Voltage (mV)')
formatForLee(gcf)
%%
for art=1%:numel(outputData.artifactData)
% art=2;
    lw = 1;
    chan =48;
    waves = [1:2:10];
    x = 0:size(outputData.artifactData(art).artifact(chan,waves,:),3)-1;
    plotX = x/30000*1000;
    figure()
%     ax1=subplot(2,1,1)
    plot(plotX,squeeze(outputData.artifactData(art).artifact(chan,waves,:)*8/1000)','linewidth',lw)
    xlabel('Time (ms)')
    ylabel('Voltage (mV)')
    title('Raw')
%     ax2=subplot(2,1,2)
%     plot(plotX,squeeze(outputDataFilteredTemp{1,1}.artifactData(art).artifact(chan,waves,:)*8/1000)','linewidth',lw)
%     xlabel('Time (ms)')
%     ylabel('Voltage (mV)')
%     linkaxes([ax1,ax2],'x')
%     title('Filtered')
% %     pause;
end
%%

%%
figDir = 'D:\Lab\Data\StimArtifact\Summary_Figures\Han_20170516_checkingChannels\';
folderpath='D:\Lab\Data\StimArtifact\Han\longFiles_20170516\';
load(strcat(folderpath,'Input_Data\','Input_structure.mat'));
load(strcat(folderpath,'Output_Data\','artifactData.mat'));
load(strcat(folderpath,'Output_Data\','chList.mat'));
load(strcat(folderpath,'Output_Data\','eList.mat'));
load(strcat(folderpath,'Output_Data\','posList.mat'));
load('Neuron_data_canonical.mat');
outputData.artifactData = artifactData; clear artifactData;
outputData.chList = chList; clear chList;
outputData.eList = eList; clear eList;
outputData.posList = posList; clear posList;

%%
lw = 1.5;
art=1;
f=figure;
f.Position = [1 41 1920 963];
for chan = 5:5:96
    subplot(4,5,chan/5)
    waves = [1:15:75];
    x = 0:size(outputData.artifactData(art).artifact(chan,waves,:),3)-1;
    plotX = x/30000*1000;
    plot(plotX,squeeze(outputData.artifactData(art).artifact(chan,waves,:)/1000)','linewidth',lw)
    ylim([-0.5 0.5])
%     xlabel('Time (ms)')
%     ylabel('Voltage (mV)')
end
chan = outputData.artifactData(art).stimChannel;
subplot(4,5,20)
plotX = x/30000*1000;
plot(plotX,squeeze(outputData.artifactData(art).artifact(chan,waves,:)*8/1000)','linewidth',lw)
fname = strcat('StimChan',num2str(chan),'_longerFile');
saveFigure(gcf,figDir,fname)

% art = 1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,19,20,21,22,23,24,25; NS = good, S = meh
% art = 11,18, NS = bad, S = bad
% art = 6,7,9,14,20,22,23 has neuron :D
% good chans = 10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,3,4,5,6,7,8,9

%% mapping: (art,stimChan) 
art = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,19,20,21,22,23,24,25];
for i = 1:numel(art)
    goodChans(i) = outputData.artifactData(art(i)).stimChannel;
end

