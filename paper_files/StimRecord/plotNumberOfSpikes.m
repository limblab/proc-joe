chan = 1;

numSpikesPost = {};
numSpikesPre = {};
for wave = 1:5
    for arr = 1:numel(arrayData)
        numSpikesPost{wave,arr} = zeros(arrayData{1}.numStims(wave),1);
        numSpikesPre{wave,arr} = zeros(arrayData{1}.numStims(wave),1);
        for st = 1:arrayData{1}.numStims(chan,wave)
            spikeIdx = find(arrayData{arr}.stimData{chan,wave} == st);
            numSpikesPost{wave,arr}(st) = numSpikesPost{wave,arr}(st) + sum(arrayData{arr}.spikeTrialTimes{chan,wave}(spikeIdx)>1.5/1000 & ...
                arrayData{arr}.spikeTrialTimes{chan,wave}(spikeIdx)<11/1000);
            numSpikesPre{wave,arr}(st) = numSpikesPre{wave,arr}(st) + sum(arrayData{arr}.spikeTrialTimes{chan,wave}(spikeIdx)<-1.5/1000 & ...
                arrayData{arr}.spikeTrialTimes{chan,wave}(spikeIdx)>-11/1000);
        end
    end
end

%%
for arrIdx = 1:size(numSpikesPost,2)
    bE = -0.5:1:100.5;
    bCPost = [];
    bCPre = [];
    figure
    hold on
    colors = {'k','r','b',[0,0.5,0],'m'};
    for wave = 1:1:5
        bCPost(wave,:) = histcounts(numSpikesPost{wave,arrIdx},bE);
        bCPost(wave,:) = bCPost(wave,:)/arrayData{1}.numStims(chan,wave);
        bCPre(wave,:) = histcounts(numSpikesPre{wave,arrIdx},bE);
        bCPre(wave,:) = bCPre(wave,:)/arrayData{1}.numStims(chan,wave);
        plot(bE(1:end-1) + mode(diff(bE))/2,bCPost(wave,:)','linewidth',1.5,'color',colors{wave})

    end
    plot(bE(1:end-1) + mode(diff(bE))/2,bCPre(wave,:)','--','linewidth',1.5,'color','k')
    xlim([0,14])
    formatForLee(gcf)
    xlabel('Number of spikes')
    ylabel('Normalized trial count')
    set(gca,'fontsize',16)

    l = legend('20\muA','30\muA','40\muA','50\muA','60\muA','baseline');
    set(l,'box','off')
end

%% do a ks test on the distribution of # spikes?

for wave = 1:5
    for arrIdx = 1:32
        [h,p(wave,arrIdx),ks2stat] = kstest2(numSpikesPre{wave,arrIdx},numSpikesPost{wave,arrIdx});
    end
end


%% heatmap?
colors = flip([0+(0:630)'/(630) 0*(0:630)' 0*(0:630)'],1);
for wave = 1:5
    figure()
    plot([1,11,11,1,1],[1,1,11,11,1],'-k','linewidth',1.5)       
    ax = gca;
    c=colormap([1,1,1]);
    ax.YTickLabel = {};
    ax.XTickLabel = {};
    for arrIdx = 1:32
        colorIdx = floor(size(colors,1)*p(wave,arrIdx));
        colorIdx = max(1,min(size(colors,1),colorIdx));
        colorToPlot = colors(colorIdx,:);
        rectangle('Position',[arrayData{arrIdx}.COL,arrayData{arrIdx}.ROW,1,1],'EdgeColor','k',...
                'FaceColor',colorToPlot,'linewidth',0.1);
    end
end

figure();
b=colorbar;
colormap(colors);
set(gca,'Visible','off');