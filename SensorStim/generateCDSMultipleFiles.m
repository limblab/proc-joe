function [] = generateCDSMultipleFiles(filepath, fileprefix)
% generates cds files for the given filepath and fileprefix
funcFolder = pwd;
files = dir([filepath fileprefix '*.mat']);
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
    end
    i=i+1;
end

% check for file existance, and if cds already exists
for i = 1:length(files)
    disp(i)
    if((exist([filepath files(i).name],'file') > 0 && ~exist([filepath files(i).name(1:end-4) '_cds.mat'],'file') ...
            && isempty(strfind(files(i).name,'EMGextra'))))
        labnum = 6;
        monkey = 'monkeyHan';
        ranBy = 'ranByRaeed';
        array = 'arrayLeftS1Area2';
        task = 'tasknone';
        cds = commonDataStructure();
        cds.file2cds([filepath files(i).name], labnum, task, monkey, ranBy, array)
        save([filepath files(i).name(1:end-4) '_cds'],'cds');
    end
end
cd(funcFolder);

end

