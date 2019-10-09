folderpath = 'C:\Users\jts3256\Desktop\Han_stim_data\Han_20190926_multielec_dukeBoardgen2\chan92\';

pwd = cd;

cd(folderpath);
stimInfoFileList = dir('*stimInfo*');

for i = 1:numel(stimInfoFileList)
    
    load([folderpath,stimInfoFileList(i).name]);
    
    % for multielectrode stim experiment
    stimInfoAdj = stimInfo;
    stimInfoAdj.waveSent =  stimInfo.waveSent(1,:)*10 + stimInfo.waveSent(2,:);
    stimInfoAdj.waveSent = reshape(repmat(stimInfoAdj.waveSent,11,1),1,numel(stimInfoAdj.waveSent)*11);
    stimInfoAdj.chanSent = reshape(repmat(stimInfoAdj.chanSent,11,1),1,numel(stimInfoAdj.chanSent)*11);
    stimInfo = stimInfoAdj;
    save([folderpath,stimInfoFileList(i).name],'stimInfo');

end










