stimChanList = [1:7,9:2:31];
chanResistorList=[1.000 10.0 1.000 10.0 1.000 10.0 1.000 2.740 2.740 2.740 3.920 3.920 3.920 4.990 4.990 4.990 5.620 5.620 5.620];
xData = ((1:size(dataStruct2.artifactData(1).artifact,3))-101)/30;
for d = 1:1%numel(dataStruct2.artifactData)
    figure();
    for art = 1:size(dataStruct2.artifactData(d).artifact,1)
        subplot(5,4,art)
        stimChan = find(stimChanList == dataStruct2.artifactData(d).stimChannel);
        data = squeeze(dataStruct2.artifactData(d).artifact(stimChan,art,:));
        dataFiltered = fliplr(data');
        plot(xData,dataFiltered)
        ylim([-300,300])
    end
    
end