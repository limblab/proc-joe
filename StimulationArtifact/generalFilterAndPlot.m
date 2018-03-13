chanList = [1:7,9:2:31];
chanResistorList=[1.000 10.0 1.000 10.0 1.000 10.0 1.000 2.740 2.740 2.740 3.920 10.0 3.920 3.920 4.990 4.990 4.990 5.620 5.620 5.620];
xData = ((1:size(dataStruct2.artifactData(1).artifact,3))-101)/30;
[b,a] = butter(6,500/(30000/2),'high');

for d = 1:numel(dataStruct2.artifactData)
    figure();
 
    subplot(2,1,1)
%     chan = 4;
%     stimChan = chan;
    stimChan = find(chanList == dataStruct2.artifactData(d).stimChannel);
    data = squeeze(dataStruct2.artifactData(d).artifact(stimChan,:,:))';
    data = [data; mean(data(end-10:end,:),1).*ones(200,20)];
    dataFiltered = fliplr(filter(b,a,fliplr(data')')')';
%     dataFiltered = data;
    plot(xData(101:end),dataFiltered(101:end-200,:))
    ylim([-250,250])
    if(chanList(stimChan) == dataStruct2.artifactData(d).stimChannel)
        title(strcat('STIMchan: ',num2str(chanList(stimChan)),' res:',num2str(chanResistorList(stimChan))));
    else
        title(strcat('chan ',num2str(chanList(stimChan)), 'res: ', num2str(chanResistorList(stimChan))));
    end
    subplot(2,1,2)
    stimChan = 4;
    
    data = squeeze(dataStruct2.artifactData(d).artifact(stimChan,:,:))';
    data = [data; mean(data(end-10:end,:),1).*ones(200,20)];
    dataFiltered = fliplr(filter(b,a,fliplr(data')')')';
%     dataFiltered = data;
    plot(xData(101:end),dataFiltered(101:end-200,:))
    ylim([-250,250])
    if(chanList(stimChan) == dataStruct2.artifactData(d).stimChannel)
        title(strcat('STIMchan: ',num2str(chanList(stimChan)),' res:',num2str(chanResistorList(stimChan))));
    else
        title(strcat('chan ',num2str(chanList(stimChan)), 'res: ', num2str(chanResistorList(stimChan))));
    end
%     suptitle(strcat('StimChan: ',num2str(dataStruct2.artifactData(d).stimChannel)))
    saveFiguresLIB(gcf,folderpath,strcat('ArtMonkey_20180128_stimChan',num2str(dataStruct2.artifactData(d).stimChannel),'_2chans_filtered'))

end