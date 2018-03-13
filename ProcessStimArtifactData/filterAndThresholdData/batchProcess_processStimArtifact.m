%% process stimulation artifacts:
pwd = cd;
folderpaths= {'R:\data\Mihili_12A3\stimRecord\Mihili_20170802_stimRecord\',...
    'R:\data\Mihili_12A3\stimRecord\Mihili_20170803_stimRecord\',...
    'R:\data\Mihili_12A3\stimRecord\Mihili_20170728_stimRecord\',...
    'R:\data\Mihili_12A3\stimRecord\Mihili_20170721_stimRecord\'};
% inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left M1 SN 1025-001452.cmp';
inputData.task='taskRW';
inputData.ranBy='ranByJoseph'; 
inputData.badChList=0;
inputData.interpulse=.000053;%in s
inputData.pWidth1=.0002;
inputData.pWidth2=.0002;
% functionName='processStimArtifactData';

MERGE_FILES = 0;
inputData.noSyncIntended = 0;
inputData.templateSubtract = 0;
inputData.templateSize = 9/1000;
inputData.blankPeriod = floor(0.3*30);
inputData.artifactDataTime = 9; % in ms

inputData.moreThanOnePulsePerWave = 0;
inputData.numPulses = 10;
inputData.pulseFrequency = 100;

inputData.thresholdMult = 3.5;

for ficker = folderpaths
    folderpath = ficker{1};
    if(~isempty(strfind(folderpath,'Han')))
        inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
        inputData.monkey = 'monkeyHan';
        inputData.array1='arrayLeftS1'; 

    elseif(~isempty(strfind(folderpath,'Mihili')))
        inputData.mapFile='mapFileR:\limblab\lab_folder\Animal-Miscellany\Mihili 12A3\Mihili Left PMd SN 6251-001460.cmp';
        inputData.monkey = 'monkeyMihili';
        inputData.array1='arrayLeftPMD'; 
    end
    disp(inputData.mapFile);
    batchList = [];
    % generates _cds and _nevData files, also writes nev file
    cd(folderpath)
    fileList = dir('*.nev');
    if(MERGE_FILES)
        endIndex = 1;
    else
        endIndex = numel(fileList);
    end

    for f = 1:endIndex
        warning('off')
%         inputData.stimsPerBump = 1;
        if(~MERGE_FILES)
            inputData.fileListExtension = fileList(f).name;
        end


        inputData.windowSize=30*10;%in points
        inputData.presample=1;%in points
        inputData.plotRange=0.5;%in mV
        inputData.lab=6;
        inputData.useSyncLabel=[];
    % dataStruct2 = runDataProcessing(functionName,folderpath,inputData)
        try
            processStimArtifactData(folderpath,inputData);
            batchList(end+1,1) = f;
        catch
            disp('error somewhere');
        end
    end
    
    try
        disp('writing nev file')
        if(~MERGE_FILES) %write nev file
            nevDataAll = [];
            totalDuration = 0;
            cd(folderpath);
            fileListNEV=dir('*_nevData*');
            fileListCDS=dir('*_cds*');
            for i = 1:numel(batchList)
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
    catch
        disp('did not write nev file, failure somewhere')
    end
end

cd(pwd);
disp('DONE -- CAN CONTINUE')
warning('on')
