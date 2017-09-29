%% set file names 

folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\dukeBoardTesting_20170929_Han\';
% folderpath = 'D:\Lab\Data\StimArtifact\Han\bumpstim\20170614\';
% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left PMd SN 6251-001460.cmp';
% mapFileName = 'C:\Users\Joseph\Desktop\Mihili Left PMd SN 6251-001460.cmp';

pwd=cd;

%% load file and parse for stim electrode number
% 1 and 4
cd(folderpath)
fileList = dir('*.mat');
fileNumber = 1;
load(fileList(fileNumber).name);
cd(pwd);

%% plot artifact data
artifactDataIdx = [19];
fields = {'breakoutData','bufferData','dukeData'};
% for field = fields
    clc;
    for idx = artifactDataIdx
        figure();
%         suptitle(field)
        for p = 25
%             subplot(2,1,p-24)
            plot(squeeze(artifactData{idx}.artifact(71:2:81,p,1:end))','r')
            hold on
            plot(squeeze(artifactData{idx}.artifact(70:2:80,p,1:end))','b')
            xlim([50,1000])
            title(num2str(p))
            formatForLee(gcf)
        end
        disp(artifactData{idx}.fileName)
%         saveFiguresLIB(gcf,folderpath,strcat(artifactData{idx}.fileName));
%         close all
    end
% end
    