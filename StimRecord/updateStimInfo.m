    folderpath = 'E:\Data\Joseph\Han_stim_data\Han_20191030_longTrains_dukeGen2\chan14\20Hz\';

    pwd = cd;

    cd(folderpath);
    stimInfoFileList = dir('*stimInfo*');


    for i = 1:numel(stimInfoFileList)
    
    load([folderpath,stimInfoFileList(i).name]);
    
    % for multielectrode stim experiment
%     stimInfoAdj = stimInfo;
 %     stimInfoAdj.waveSent =  stimInfo.waveSent(1,:)*10 + stimInfo.waveSent(2,:);
%     stimInfoAdj.waveSent = reshape(repmat(stimInfoAdj.waveSent,11,1),1,numel(stimInfoAdj.waveSent)*11);
%     stimInfoAdj.chanSent = reshape(repmat(stimInfoAdj.chanSent,11,1),1,numel(stimInfoAdj.chanSent)*11);
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


%% fix stim info if someone (Juliet) messed up

folderpath = 'C:\Users\joh8881\Desktop\Han_20190924_trains_noAmp\Chan66\';

pwd = cd;

cd(folderpath);
outputDataFileList = dir('*outputData*');
stimInfoFileList = dir('*stimInfo*');

for i = 1:numel(outputDataFileList)
    
    load([folderpath,outputDataFileList(i).name]);
    
    stimInfo = outputData.stimInfo;
    stimInfo.chanSent = outputData.waveforms.chanSent';
    stimInfo.waveSent = outputData.waveforms.waveSent';
    stimInfo.parameters = outputData.waveforms.parameters;
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






