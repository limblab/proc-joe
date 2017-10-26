%% set file names 
% 
% folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20171010_dukeBoard\';
% % folderpath = 'D:\Lab\Data\StimArtifact\Han\bumpstim\20170614\';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
% % mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left PMd SN 6251-001460.cmp';
% % mapFileName = 'C:\Users\Joseph\Desktop\Mihili Left PMd SN 6251-001460.cmp';
% 
% pwd=cd;
% 
% %% load file and parse for stim electrode number
% % 1 and 4
% cd(folderpath)
% fileList = dir('*.mat');
% fileNumber = 1;
% load(fileList(fileNumber).name);
% cd(pwd);

%% plot artifact data
artifactDataIdx = [1];
fields = {'breakoutData','bufferData','dukeData'};
% for field = fields
    clc;
    for idx = artifactDataIdx
        figure();
%         suptitle(field)
        for p = 73
%             subplot(2,1,p-24)
            xData = (0:(size(artifactData(idx).artifact,3)-1))/30 - 3.;
            plot(xData,squeeze(artifactData(idx).artifact(p,1:10:end,1:end))','r')
            hold on
            plot(xData,squeeze(artifactData(idx).artifact(p,2:10:end,1:end))','b')
            xlim([-0.3,8])
%             title(num2str(p))
            formatForLee(gcf)

        end
%         disp(artifactData(idx).fileName)
        saveFiguresLIB(gcf,folderpath,fileName);
%         close all
    end
% end
    