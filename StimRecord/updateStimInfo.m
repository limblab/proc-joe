folderpath = 'C:\Users\jts3256\Desktop\Han_stim_data\Han_20190923_trains_noAmp\chan21\';

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










