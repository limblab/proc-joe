chan = 1;

numSpikesPost = {};
numSpikesPre = {};
for wave = 1:5
    numSpikesPost{wave} = zeros(arrayData{1}.numStims(chan,wave),1);
    numSpikesPre{wave} = zeros(arrayData{1}.numStims(chan,wave),1);
    for st = 1:arrayData{1}.numStims(chan,wave)
        for arr = 1:numel(arrayData)
            spikeIdx = find(arrayData{arr}.stimData{chan,wave} == st);
            numSpikesPost{wave}(st) = numSpikesPost{wave}(st) + sum(arrayData{arr}.spikeTrialTimes{chan,wave}(spikeIdx)>1/1000 & ...
                arrayData{arr}.spikeTrialTimes{chan,wave}(spikeIdx)<11/1000);
            numSpikesPre{wave}(st) = numSpikesPre{wave}(st) + sum(arrayData{arr}.spikeTrialTimes{chan,wave}(spikeIdx)<-1/1000 & ...
                arrayData{arr}.spikeTrialTimes{chan,wave}(spikeIdx)>-11/1000);
        end
    end
end

%%
bE = -0.5:1:100.5;
bCPost = [];
bCPre = [];
figure
hold on
colors = {'k','r','b',[0,0.5,0],'m'};
for wave = 1:1:5
    bCPost(wave,:) = histcounts(numSpikesPost{wave},bE);
    bCPost(wave,:) = bCPost(wave,:)/arrayData{1}.numStims(chan,wave);
    bCPre(wave,:) = histcounts(numSpikesPre{wave},bE);
    bCPre(wave,:) = bCPre(wave,:)/arrayData{1}.numStims(chan,wave);
    plot(bE(1:end-1) + mode(diff(bE))/2,bCPost(wave,:)','linewidth',1.5,'color',colors{wave})

end
plot(bE(1:end-1) + mode(diff(bE))/2,bCPre(wave,:)','--','linewidth',1.5,'color','k')
xlim([0,14])
formatForLee(gcf)
xlabel('Number of spikes')
ylabel('Normalized trial count')
set(gca,'fontsize',16)

l = legend('20\muA','30\muA','40\muA','50\muA','60\muA','baseline')
set(l,'box','off')