folderpath = 'C:\Users\joh8881\Desktop\Han_20190924_trains_noAmp\Chan1\';

pwd = cd;

cd(folderpath);
stimInfoFileList = dir('*stimInfo*');

for i = 1:numel(stimInfoFileList)
    
    load([folderpath,stimInfoFileList(i).name]);
%     stimInfoAdj = stimInfo;
%     stimInfoAdj.waveSent = reshape(repmat(stimInfo.waveSent,11,1),1,numel(stimInfo.waveSent)*11);
%     stimInfoAdj.chanSent = reshape(repmat(stimInfo.chanSent,11,1),1,numel(stimInfo.chanSent)*11);
%     stimInfo = stimInfoAdj;
%     save([folderpath,stimInfoFileList(i).name],'stimInfo');

end


%% undo

folderpath = 'C:\Users\joh8881\Desktop\Han_20190924_trains_noAmp\Chan66\';

pwd = cd;

cd(folderpath);
stimInfoFileList = dir('*stimInfo*');

for i = 1:numel(stimInfoFileList)
    
    load([folderpath,stimInfoFileList(i).name]);
    stimInfoTemp = stimInfo;
    stimInfoUndo.waveSent = [];
    stimInfoUndo.chanSent = [];
    index = numel(stimInfoTemp.waveSent)
    for n = 1:1100
        stimInfoUndo.waveSent = [stimInfoUndo.waveSent,stimInfoTemp.waveSent(n*11)];
        stimInfoUndo.chanSent = [stimInfoUndo.chanSent,stimInfoTemp.chanSent(n*11)];
    end
    stimInfo = stimInfoUndo;
    save([folderpath,stimInfoFileList(i).name],'stimInfo');

end

%% undo joe

folderpath = 'C:\Users\joh8881\Desktop\Han_20190924_trains_noAmp\Chan1\';

pwd = cd;

cd(folderpath);
stimInfoFileList = dir('*stimInfo*');

for i = 1:numel(stimInfoFileList)
    
    load([folderpath,stimInfoFileList(i).name]);
    stimInfo.waveSent = stimInfo.waveSent(1:11:end);
    stimInfo.chanSent = stimInfo.chanSent(1:11:end);
    save([folderpath,stimInfoFileList(i).name],'stimInfo');

end





