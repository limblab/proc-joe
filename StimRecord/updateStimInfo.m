    folderpath = 'E:\Data\Joseph\Han_stim_data\Han_20191030_longTrains_dukeGen2\chan14\20Hz\';

    pwd = cd;

    cd(folderpath);
    stimInfoFileList = dir('*stimInfo*');

    for i = 1:numel(stimInfoFileList)

        load([folderpath,stimInfoFileList(i).name]);

        stimInfoAdj = stimInfo;
        % for multielectrode stim experiment
    %     stimInfoAdj.waveSent =  stimInfo.waveSent(1,:)*10 + stimInfo.waveSent(2,:);
        num_rep = numel(stimInfoAdj.stimOn)/numel(stimInfoAdj.chanSent);
        
        stimInfoAdj.waveSent = reshape(repmat(stimInfoAdj.waveSent,num_rep,1),1,numel(stimInfoAdj.waveSent)*num_rep);
        stimInfoAdj.chanSent = reshape(repmat(stimInfoAdj.chanSent,num_rep,1),1,numel(stimInfoAdj.chanSent)*num_rep);
        stimInfo = stimInfoAdj;
        save([folderpath,stimInfoFileList(i).name],'stimInfo');

    end



%% undo change

    folderpath = 'E:\Data\Joseph\Han_stim_data\Han_20191029_longTrains_dukeGen2\chan21\';

    pwd = cd;

    cd(folderpath);
    stimInfoFileList = dir('*stimInfo*');

    for i = 1:numel(stimInfoFileList)

        load([folderpath,stimInfoFileList(i).name]);

        stimInfoAdj = stimInfo;
        % for multielectrode stim experiment
    %     stimInfoAdj.waveSent =  stimInfo.waveSent(1,:)*10 + stimInfo.waveSent(2,:);

    %     stimInfoAdj.waveSent = stimInfoAdj.waveSent(1:11:end);
    %     stimInfoAdj.chanSent = stimInfoAdj.chanSent{1:11:end};
        stimInfoAdj.chanSent = mat2cell(stimInfoAdj.chanSent,1,ones(1,numel(stimInfoAdj.chanSent)));
    %     stimInfo = stimInfoAdj;
        save([folderpath,stimInfoFileList(i).name],'stimInfo');

    end






