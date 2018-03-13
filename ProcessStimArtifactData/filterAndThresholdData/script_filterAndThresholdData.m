clear%% process stimulation artifacts:
pwd = cd;
folderpath= 'R:\data\Han_13B1\Raw\Han_20180304_chic201802\normalSetup\';
% inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left PMd SN 6251-001460.cmp'; % chips mapfile location
inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
% inputData.mapFile = 'mapFileR:\limblab\lab_folder\Animal-Miscellany\Chips_12H1\map_files\left S1\SN 6251-001455.cmp';
inputData.task='taskCO';
inputData.ranBy='ranByJoseph'; 
inputData.array1='arrayLeftS1'; 
inputData.monkey='monkeyChips';

inputData.dukeBoardChannel = -1;
inputData.dukeBoardLabel = 'ainp15';

inputData.issueExists = 0;

inputData.badChList=0;
inputData.interpulse=.000053;%in s
inputData.pWidth1=.0002;
inputData.pWidth2=.0002;
% functionName='processStimArtifactData';
MERGE_FILES = 0;
inputData.noSyncIntended = 0;
inputData.templateSubtract = 0;
inputData.templateSize = 99/1000;

inputData.blankPeriod = floor(0.0*30);
inputData.artifactDataTime = 10; % in ms

inputData.preOffset = 22;
inputData.postOffset = 25;

inputData.moreThanOnePulsePerWave = 0;
inputData.numPulses = 10;
inputData.pulseFrequency = 100;


inputData.thresholdMult = 3.5;
inputData.artifactSkip = 1;

inputData.maxChunkLength = 5000*30; % 5 second chunk maximum
%% generates _cds and _nevData files, also writes nev file
cd(folderpath)
fileList = dirSorted('*.nev');

if(MERGE_FILES)
    endIndex = 1;
else
    endIndex = numel(fileList);
end

% save input data
save(strcat(fileList(1).name(1:end-4),'_inputData.mat'),'inputData');

artifactData = cell(endIndex,1);
% process data
for f = 3:endIndex
    warning('off')
    inputData.stimsPerBump = 1;
    if(~MERGE_FILES)
        inputData.fileListExtension = fileList(f).name;
    end
    
    inputData.windowSize=30*10;%in points
    inputData.presample=1;%in points
    inputData.plotRange=0.5;%in mV
    inputData.lab=6;
    inputData.useSyncLabel=[];
% dataStruct2 = runDataProcessing(functionName,folderpath,inputData)
    processStimArtifactData(folderpath,inputData);
%     artifactData{f}.artifact = extractArtifactData(folderpath,inputData);
%     artifactData{f}.fileName = fileList(f).name;
end

% save('','artifactData')

% write nev file
% disp('writing nev file')
% if(~MERGE_FILES) %write nev file
%     nevDataAll = [];
%     totalDuration = 0;
%     cd(folderpath);
%     fileListNEV=dirSorted('*_nevData*');
%     fileListCDS=dirSorted('*_cds.mat*');
%     for i = 1:numel(fileListNEV)
%         load(fileListNEV(i).name);
%         load(fileListCDS(i).name);
%         if(i==1)
%             nevDataAll.ts = nevData.ts;
%             nevDataAll.waveforms = nevData.waveforms(:,:);
%             nevDataAll.elec = nevData.elec(:,:);
%         else
%             nevDataAll.ts(end+1:end+numel(nevData.ts),:) = nevData.ts + totalDuration;
%             nevDataAll.waveforms(end+1:end+numel(nevData.ts),:) = nevData.waveforms(:,:);
%             nevDataAll.elec(end+1:end+numel(nevData.ts),:) = nevData.elec(:,:);
%         end
%         disp(num2str(cds.meta.duration - max(nevData.ts)))
%         
%         totalDuration = totalDuration + cds.meta.duration;
%     end
%     packetWidth = 104;
%     filename = strcat(fileListNEV(1).name(1:12),'_merged');
%     mapFilename = inputData.mapFile(8:end);
%     comments = '';
%     writeNEV(nevDataAll, packetWidth, filename, mapFilename, comments )
%     cd(pwd);
% else
%     fileListNEV=dirSorted('*_nevData*');
%     load(fileListNEV(1).name);
%     packetWidth = 104;
%     filename = strcat(fileListNEV(1).name(1:12),'_merged');
%     mapFilename = inputData.mapFile(8:end);
%     comments = '';
%     writeNEV(nevData, packetWidth, filename, mapFilename, comments )
%     cd(pwd);
% end

cd(pwd);
disp('DONE -- CAN CONTINUE')
warning('on')

%% writeNev from _nevData file if failed above
% inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left PMd SN 6251-001460.cmp'; % chips mapfile location
nevDataAll = [];
cdsAll = [];
totalDuration = 0;
pwd = cd;
cd(folderpath);
fileListNEV=dirSorted('*_nevData*');
fileListCDS = dirSorted('*_cds.mat*');
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
    disp(num2str(cds.meta.duration - max(nevData.ts)))
    totalDuration = totalDuration + cds.meta.duration;
end
packetWidth = 104;
filename = strcat(fileListNEV(1).name(1:14),'_merged');
mapFilename = inputData.mapFile(8:end);
comments = '';
writeNEV(nevDataAll, packetWidth, filename, mapFilename, comments )
cd(pwd);
disp('done, can continue')
%% sort *_merged and call it *_merged-s

%% load in *_merged-s, save the units and move that to cds 
pwd = cd;
cd(folderpath);
disp('started')
if(MERGE_FILES)
    fileListCDS = dirSorted('*_cds.mat*');
    load(fileListCDS(1).name);
    fileList = dirSorted('*_merged-s*');
    filename = fileList(1).name;
    labnum = 6;
    monkey = inputData.monkey;
    ranBy = inputData.ranBy;
    array = inputData.array1;
    task = inputData.task;
    spikes = commonDataStructure();
    spikes.file2cds([folderpath,fileList(1).name],ranBy,array,monkey,labnum,task,inputData.mapFile);
    cds.units = spikes.units;
    clear spikes
    cd(folderpath)
    save(strcat(filename(1:end-13),'_processed'),'cds','-v7.3');
    cd(pwd);
elseif(~MERGE_FILES)
    fileListCDS = dirSorted('*_cds.mat*');
    nevFileList = dirSorted('*_merged-s.NEV*');
    nevFile = nevFileList(1).name;
    labnum = 6;
    monkey = inputData.monkey;
    ranBy = inputData.ranBy;
    array = inputData.array1;
    task = inputData.task;
    spikes = commonDataStructure();
    spikes.file2cds([folderpath,nevFileList(1).name],ranBy,array,monkey,labnum,'ignoreJumps',task,inputData.mapFile);
    
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

disp('done with this step')

%% if need to merge cds's now, do so here
% folderpath = 'R:\data\Mihili_12A3\stimRecord\Mihili_20170717_stimRecord\';
pwd = cd;
cd(folderpath)
fileListProcessed = dirSorted('*_cds_processed.mat');
cdsAll = [];

for f = 1:numel(fileListProcessed)
    load(fileListProcessed(f).name);
    % trim artifact data to remove excess entries
    lastArtifactIdx = -1;
    for artIdx = 1:size(cds.artifactData.artifact,1)
        if(lastArtifactIdx == -1 && sum(sum(squeeze(cds.artifactData.artifact(artIdx,:,:)) == zeros(size(cds.artifactData.artifact,2),size(cds.artifactData.artifact,3)))) == size(cds.artifactData.artifact,2)*size(cds.artifactData.artifact,3))
            lastArtifactIdx = artIdx;
            disp(lastArtifactIdx)
        end
    end
    cds.artifactData.artifact = cds.artifactData.artifact(1:lastArtifactIdx-1,:,:);
    cds.artifactData.t = cds.artifactData.t(1:lastArtifactIdx-1);
    
    % merge files into 1 cds
    if(f==1)
        % rewrite everything into cds -- which is really not a
        % commonDataStructure object but whatever
        cdsAll.offsets = [0];
        cdsAll.meta = cds.meta;
        cdsAll.meta.hasLfp = 0;
        if(isfield(cds,'kin'))
            cdsAll.kin = cds.kin;
        else
            cdsAll.kin = [];
        end
        cdsAll.force = cds.force;
        cdsAll.lfp = {};
        cdsAll.emg = cds.emg;
        cdsAll.analog = {};
        cdsAll.stimOn = cds.stimOn;
        cdsAll.stimOff = cds.stimOff;
        cdsAll.triggers = cds.triggers;
        cdsAll.units = cds.units;
        if(any(isfield(cds,'trials')))
            cdsAll.trials = cds.trials;
        end
        cdsAll.aliasList = cds.aliasList;
        cdsAll.operationLog = cds.operationLog;
        cdsAll.meta.duration = cds.meta.duration;
        cdsAll.artifactData.t = cds.artifactData.t(:);
        cdsAll.artifactData.artifact = cds.artifactData.artifact(:,:,:);
        cdsAll.rawData = cds.rawData;
    else     
        % port over kin data
        if(isfield(cds,'kin') && ~isempty(cds.kin))
            tempKin = cds.kin;
            tempKin.t = tempKin.t + cdsAll.meta.duration;
            cdsKin = cdsAll.kin;
            cdsKin{end+1:end+size(tempKin,1),:} = tempKin{:,:};
            cdsAll.kin = cdsKin;
        end
        % port over artifact data
        cdsAll.artifactData.t(size(cdsAll.artifactData.artifact,1)+1:size(cdsAll.artifactData.artifact,1)+size(cds.artifactData.artifact,1)) = cds.artifactData.t(:) + cdsAll.meta.duration;
        cdsAll.artifactData.artifact(size(cdsAll.artifactData.artifact,1)+1:size(cdsAll.artifactData.artifact,1)+size(cds.artifactData.artifact,1),:,:) = cds.artifactData.artifact(:,:,1:size(cdsAll.artifactData.artifact,3));
        % port over stimOn and stimOff data
        cdsAll.stimOn = [cdsAll.stimOn;cds.stimOn + cdsAll.meta.duration];
        cdsAll.stimOff = [cdsAll.stimOff;cds.stimOff + cdsAll.meta.duration];
        % port over raw data
        cdsAll.rawData.ts(end+1:end+numel(cds.rawData.ts)) = cds.rawData.ts + cdsAll.meta.duration;
        cdsAll.rawData.elec(end+1:end+numel(cds.rawData.ts)) = cds.rawData.elec;
        cdsAll.rawData.waveforms(end+1:end+numel(cds.rawData.ts),:) = cds.rawData.waveforms;
        % port over trial data
        if(any(isfield(cdsAll,'trials')))
            try
                tempTrials = cds.trials;
                tempTrials.number = tempTrials.number + cdsAll.trials.number(end);
                tempTrials.startTime = tempTrials.startTime + cdsAll.meta.duration;
                tempTrials.endTime = tempTrials.endTime + cdsAll.meta.duration;
                tempTrials.tgtOnTime = tempTrials.tgtOnTime + cdsAll.meta.duration;
                tempTrials.goCueTime = tempTrials.goCueTime + cdsAll.meta.duration;
                tempTrials.bumpTime = tempTrials.bumpTime + cdsAll.meta.duration;
                tempTrials.stimTime = tempTrials.stimTime + cdsAll.meta.duration;
                cdsTrials = cdsAll.trials;
                cdsTrials(end+1:end+size(tempTrials,1),:) = tempTrials(:,:);
                cdsAll.trials = cdsTrials;
                cdsAll.meta.numTrials = cdsAll.meta.numTrials + cds.meta.numTrials;
                cdsAll.meta.numReward = cdsAll.meta.numReward + cds.meta.numReward;
                cdsAll.meta.numAbort = cdsAll.meta.numAbort + cds.meta.numAbort;
                cdsAll.meta.numFail = cdsAll.meta.numFail + cds.meta.numFail;
                cdsAll.meta.numIncomplete = cdsAll.meta.numIncomplete + cds.meta.numIncomplete;
            catch
            end
        end
        
        % update units information
        for nn = 1:size(cds.units,2)
            nnAll = 1;
            cds.units(nn).spikes{:,1} = cds.units(nn).spikes{:,1} + cdsAll.meta.duration;
            while(nnAll <= size(cds.units,2) && (cds.units(nnAll).chan ~= cds.units(nn).chan || cds.units(nnAll).ID ~= cds.units(nn).ID))
                nnAll = nnAll + 1;
            end
            if(nnAll <= size(cds.units,2))
                cdsAll.units(nnAll).spikes{end+1:end+size(cds.units(nn).spikes,1),:} = cds.units(nn).spikes{:,:};
            end
            
        end
        
        % update meta information
        cdsAll.offsets(end+1,1) = cdsAll.meta.duration;
        cdsAll.meta.dataWindow(2) = cds.meta.dataWindow(2) + cds.meta.dataWindow(2);
        cdsAll.meta.duration = cdsAll.meta.duration + cds.meta.duration;
    end
    
end

cds = cdsAll;
clear cdsAll;
save(strcat(fileListProcessed(1).name(1:26),'_all_processed'),'cds','-v7.3');

cd(pwd)
disp('done merging')
%% if an interleaved set of trials, merge waveform sent information into cds
cd(folderpath);
fileListCDS = dirSorted('*_processed.mat*');
% load(fileListCDS(1).name);
cds.waveforms.waveSent = [];

fileListWaves = dirSorted('*_waveformsSent_*');
for f = 1:numel(fileListWaves)
    load(fileListWaves(f).name);
    load(fileListCDS(f).name)
    disp(num2str(f))
    d=find(diff(cds.stimOn)>10);
    if(~isempty(d))
        if(f==1)
            cds.waveforms.waveSent(1:d(1),1) = waveforms.waveSent(1:d(1));
            cds.waveforms.chanSent(1:d(1),1) = waveforms.chanSent(1:d(1));
        elseif(f==numel(fileListWaves ))
            numWaves = numel(cds.stimOn) - numel(cds.waveforms.waveSent);
            cds.waveforms.waveSent(end+1:end+numWaves) = waveforms.waveSent(1:numWaves);
            cds.waveforms.chanSent(end+1:end+numWaves) = waveforms.chanSent(1:numWaves);
        else
            numWaves = d(f)-d(f-1);
            cds.waveforms.waveSent(end+1:end+numWaves) = waveforms.waveSent(1:numWaves);
            cds.waveforms.chanSent(end+1:end+numWaves) = waveforms.chanSent(1:numWaves);
        end
    else
        cds.waveforms.waveSent = waveforms.waveSent;
    end
    cds.waveforms.parameters = waveforms.parameters;
    save(fileListCDS(f).name,'cds','-v7.3');
end

% save(fileListCDS(1).name,'cds','-v7.3');

disp('done interleaving')