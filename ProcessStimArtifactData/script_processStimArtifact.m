%% process stimulation artifacts:
pwd = cd;
folderpath= 'R:\data\Han_13B1\Raw\20170629\chan42stim\';
% functionName='processStimArtifactData';
MERGE_FILES = 1;
%% generates _cds and _nevData files
cd(folderpath)
fileList = dir('*.nev');
if(MERGE_FILES)
    endIndex = 1;
else
    endIndex = numel(fileList);
end

for f = 1:endIndex
    warning('off')
    inputData.stimsPerBump = 1;
    if(~MERGE_FILES)
        inputData.fileListExtension = fileList(f).name;
    end
    inputData.task='taskCObump';
    inputData.ranBy='ranByJoseph'; 
    inputData.array1='arrayS1'; 
    inputData.monkey='monkeyHan';
    inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp'; % chips mapfile location
    inputData.badChList=0;
    inputData.interpulse=.000053;%in s
    inputData.pWidth1=.0002;
    inputData.pWidth2=.0002;

    inputData.windowSize=30*10;%in points
    inputData.presample=1;%in points
    inputData.plotRange=0.5;%in mV
    inputData.lab=6;
    inputData.useSyncLabel=[];

% dataStruct2 = runDataProcessing(functionName,folderpath,inputData)
    processStimArtifactData(folderpath,inputData);
    
end

disp('writing nev file')
if(~MERGE_FILES) %write nev file
    nevDataAll = [];
    totalDuration = 0;
    cd(folderpath);
    fileListNEV=dir('*_nevData*');
    fileListCDS=dir('*_cds*');
    for i = 1:numel(fileList)
        load(fileListNEV(i).name);
        load(fileListCDS(i).name);
        if(i==1)
            nevDataAll.ts = nevData.ts;
            nevDataAll.waveforms = nevData.waveforms(:,:);
            nevDataAll.elec = nevData.elec(:,:);
        else
            nevDataAll.ts(end+1:end+numel(nevData.ts),:) = nevData.ts;
            nevDataAll.waveforms(end+1:end+numel(nevData.ts),:) = nevData.waveforms(:,:);
            nevDataAll.elec(end+1:end+numel(nevData.ts),:) = nevData.elec(:,:);
        end

        totalDuration = totalDuration + cds.meta.duration;
    end
    packetWidth = 104;
    filename = strcat(fileListNEV(1).name(1:12),'_merged');
    mapFilename = inputData.mapFile(8:end);
    comments = '';
    writeNEV(nevDataAll, packetWidth, filename, mapFilename, comments )
    cd(pwd);
else
    fileListNEV=dir('*_nevData*');
    load(fileListNEV(1).name);
    packetWidth = 104;
    filename = strcat(fileListNEV(1).name(1:12),'_merged');
    mapFilename = inputData.mapFile(8:end);
    comments = '';
    writeNEV(nevData, packetWidth, filename, mapFilename, comments )
    cd(pwd);
end


cd(pwd);
disp('DONE -- CAN CONTINUE')
warning('on')

%% writeNev from _nevData file if failed above

nevDataAll = [];
cdsAll = [];
totalDuration = 0;
pwd = cd;
cd(folderpath);
fileListNEV=dir('*_nevData*');
fileListCDS = dir('*_cds.mat*');
for i = 1:numel(fileListNEV)
    load(fileListNEV(i).name);
    load(fileListCDS(i).name);

    if(i==1)
        nevDataAll.ts = nevData.ts;
        nevDataAll.waveforms = nevData.waveforms(:,:);
        nevDataAll.elec = nevData.elec(:,:);
    else
        nevDataAll.ts(end+1:end+numel(nevData.ts),:) = nevData.ts + totalDuration;
        nevDataAll.waveforms(end+1:end+numel(nevData.ts),:) = nevData.waveforms(:,:);
        nevDataAll.elec(end+1:end+numel(nevData.ts),:) = nevData.elec(:,:);
    end
      
    totalDuration = totalDuration + cds.meta.duration;
end
packetWidth = 104;
filename = strcat(fileListNEV(1).name(1:12),'_merged');
mapFilename = inputData.mapFile(8:end);
comments = '';
writeNEV(nevDataAll, packetWidth, filename, mapFilename, comments )
cd(pwd);
%% sort *_merged and call it *_merged-s

%% load in *_merged-s, save the units and move that to cds 
pwd = cd;
cd(folderpath);
if(MERGE_FILES)
    fileListCDS = dir('*_cds.mat*');
    load(fileListCDS(1).name);
    fileList = dir('*_merged-s.NEV*');
    filename = fileList(1).name;
    labnum = 6;
    monkey = 'monkeyHan';
    ranBy = 'ranByJoseph';
    array = 'arrayLeftS1';
    task = 'taskRW';
    spikes = commonDataStructure();
    spikes.file2cds([folderpath,filename],ranBy,array,monkey,labnum,'ignoreJumps',task);
    cds.units = spikes.units;
    clear spikes
    cd(folderpath)
    save(strcat(filename(1:end-13),'_processed'),'cds','-v7.3');
    cd(pwd);
elseif(~MERGE_FILES)
    nevFileList = dir('*_merged-s.NEV*');
    nevFile = nevFileList(1).name;
    labnum = 6;
    monkey = 'monkeyHan';
    ranBy = 'ranByJoseph';
    array = 'arrayLeftS1';
    task = 'taskRW';
    spikes = commonDataStructure();
    spikes.file2cds([folderpath,nevFile],ranBy,array,monkey,labnum,'ignoreJumps',task);
    
    % split spikes units and store in the cds
    totalDuration = 0;
    for f = 1:numel(fileListCDS)
        filename = fileListCDS(f).name;
        load(filename);
        units = spikes.units;
        
        for nn = 1:size(spikes.units,2)
            units(nn).spikes.ts = units(nn).spikes.ts-totalDuration;
            spikeMask = units(nn).spikes.ts > 0 & ...
                units(nn).spikes.ts < cds.meta.duration;
            units(nn).spikes = units(nn).spikes(spikeMask,:);
        end
        cds.units = units;
        save(strcat(filename(1:end-4),'_processed'),'cds','-v7.3');
        totalDuration = totalDuration + cds.meta.duration;
    end
    
end

%% if an interleaved set of trials, merge waveform sent information into cds
if(MERGE_FILES)
    fileListCDS = dir('*_processed.mat*');
    load(fileListCDS(1).name);
    cds.waveforms.waveSent = [];
    
    fileListWaves = dir('*_waveformsSent_*');
    for f = 1:numel(fileListWaves)
        load(fileListWaves(f).name);
        d=find(diff(cds.stimOn)>1);
        if(f==1)
            cds.waveforms.waveSent(1:d(1),1) = waveforms.waveSent(1:d(1));
        elseif(f==numel(fileListWaves))
            numWaves = numel(cds.stimOn) - numel(cds.waveforms.waveSent);
            cds.waveforms.waveSent(end+1:end+numWaves) = waveforms.waveSent(1:numWaves);
        else
            numWaves = d(f)-d(f-1);
            cds.waveforms.waveSent(end+1:end+numWaves) = waveforms.waveSent(1:numWaves);
        end
        cds.waveforms.parameters = waveforms.parameters;
        
    end
    
    save(fileListCDS(1).name,'cds','-v7.3');
else
    
end

