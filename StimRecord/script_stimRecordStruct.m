%% this script combines information from one data session into a struct that
%% data from all data sessions

srs = stimRecordStruct();

folderpaths = {'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Mihili_20170721\',...
    'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Mihili_20170720\'};

filenames = {'Mihili_20170721_all_processed.mat',...
    'Mihili_20170720_chanINTERL_all_processed.mat'};

for i = 1:numel(folderpaths)
    srs.addData(folderpaths{i},filenames{i});
end


