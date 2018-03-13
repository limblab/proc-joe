%% process stimulation artifacts:
pwd = cd;
folderpath='C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20180311_stimswitch\';
functionName='processStimArtifact_filter'; % includes acausal filtered data

%%
inputData.task='taskCObump';
inputData.ranBy='ranByTucker'; 
inputData.array1='arrayS1'; 
inputData.monkey='monkeyHan';
inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp'; % mapfile location

inputData.badChList=[];
inputData.interpulse=.000053;%in s
inputData.pWidth1=.0002;
inputData.pWidth2=.0002;

inputData.windowSize=30*18;%in points
inputData.presample=100;%in points
inputData.plotRange=0.2;%in mV
inputData.lab=6;
inputData.useSyncLabel=[];

output_data = runDataProcessing(functionName,folderpath,inputData);
cd(pwd);


%% plot artifact data on the stimulated channel only
for artIdx = 1:numel(output_data.artifactData)
    figure()
    for chanIdx= 1:size(output_data.artifactData(artIdx).artifact,1)
        xData = ((1:size(output_data.artifactData(artIdx).artifact,3))-1)/30 - inputData.presample/30;
        subplot(6,6,chanIdx)
%         plot(xData,squeeze(output_data.artifactData(artIdx).artifactFiltered(output_data.artifactData(artIdx).stimChannel,1:10:50,:))','r');
%         hold on
%         plot(xData,squeeze(output_data.artifactData(artIdx).artifactFiltered(output_data.artifactData(artIdx).stimChannel,2:10:50,:))','b');
        plot(xData,squeeze(output_data.artifactData(artIdx).artifact(chanIdx,1:10:50,:))','r');
        hold on
        plot(xData,squeeze(output_data.artifactData(artIdx).artifact(chanIdx,2:10:50,:))','b');
        ylim([-10000,10000])
        xlim([0,5])
        formatForLee(gcf);
%         title(fileList(artIdx).name(17:32))
%         saveFiguresLIB(gcf,folderpath,strcat(fileList(artIdx).name(1:end-4),'_allChannels'));
    end
end

