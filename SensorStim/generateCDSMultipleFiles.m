function [] = generateCDSMultipleFiles(filepath, fileprefix)
% generates cds files for the given filepath and fileprefix
funcFolder = pwd;
filesAll = dir([filepath fileprefix '*.mat']);
files = filesAll;
% remove '_cds' files and '-s' files
i = 1;
while i <= length(files)
    if(~isempty(strfind(files(i).name,'-s')))
        files(i) = []; % removes row i from struct files
        i=i-1;
    elseif(~isempty(strfind(files(i).name,'_cds')))
        files(i) = []; % removes row i from struct files
        i=i-1;
    elseif(~isempty(strfind(files(i).name,'EMGextra')))
        files(i) = []; % removes row i from struct files
        i=i-1;
    elseif(~isempty(strfind(files(i).name,'Sweep')))
        files(i) = []; % removes row i from struct files
        i=i-1;
    elseif(~isempty(strfind(lower(files(i).name),'neurons')))
        files(i) = [];
        i=i-1;
    elseif(~isempty(strfind(lower(files(i).name),'rightcuneate')))
        files(i) = [];
        i=i-1;
    end
    i=i+1;
end

% check for file existance, and if cds already exists
for i = 1:numel(files)
    if((exist([filepath files(i).name],'file') > 0 && ~exist([filepath files(i).name(1:end-4) '_cds.mat'],'file') ...
            && isempty(strfind(files(i).name,'EMGextra'))))
        labnum = 6;
        monkey = 'monkeyLando';
        ranBy = 'ranByChris';
        array = 'arrayLeftS1';
        task = 'taskRW';
        cds = commonDataStructure();
        cds.file2cds([filepath files(i).name], labnum, task, monkey, ranBy, array)
        % check for EMGextra file
        fileExtra = dir([filepath files(i).name(1:end-17) '_EMGextra*.mat']);
        if(numel(fileExtra) > 0)
            cds.file2cds([filepath fileExtra(1).name],labnum,task,monkey,ranBy,array);
        end
        % deal with dual recordings
        % muscle name
        underscores = strfind(files(i).name,'_');
        muscleName = files(i).name(underscores(4)+1:underscores(5)-1);
        for f = 1:numel(filesAll)
            if(~isempty(strfind(lower(filesAll(f).name),'rightcuneate')) && ...
                ~isempty(strfind(filesAll(f).name,muscleName)))
                array = 'arrayRightCuneate';
                cds.file2cds([filepath filesAll(f).name],labnum,task,monkey,ranBy,array);
            end
        end
        save([filepath files(i).name(1:end-4) '_cds'],'cds');
    end
end
cd(funcFolder);

end

