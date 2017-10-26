function [] = extractSpikingInfo(folderpath,inputData)

%script to load stimulation files and generate perievent plots of 30khz
    %data. Formatted to work with runDataProcessing
    outputFigures = [];
    outputData=[];
    %get list of all files in the folder:
    templateSize = inputData.templateSize*30000;
    noSync = 0;
    noSyncIntended = inputData.noSyncIntended;
    
    if ~strcmp(folderpath(end),filesep)
        folderpath=[folderpath,filesep];
    end
    cd(folderpath);
    if(any(isfield(inputData,'fileListExtension')))
        fileList=dir(inputData.fileListExtension);
    else
        fileList=dir('*.nev');
    end
    mapData=loadMapFile(inputData.mapFile(8:end));
    for i=1:size(mapData,1)
        mapData.label{i}=[inputData.array1(6:end),mapData.label{i}];
    end
    %establish trackers for information that is common to all files:
    eList={};
    posList=[];
    chList=[];
    % build filter
    [bFilter,aFilter] = butter(6,500/(30000/2),'high');
    filterParams.b = bFilter;
    filterParams.a = aFilter;
    filterParams.order = 6;
    filterParams.type = 'high';
    thresholdMult = inputData.thresholdMult;
    % variables to store spike information   
%     preOffset = 27;
%     postOffset = 20;
    preOffset = inputData.preOffset;
    postOffset = inputData.postOffset;
    maxAmplitude = 700; %uV
    lengthWave = preOffset+postOffset+1;
    numZeros = 200;
    
    % cds and extraseconds for merge purposes
    cds = [];
    nevData = [];
    rawData = [];
    artifactDataTime = inputData.artifactDataTime; % in ms
    artifactData.artifact = zeros(3000,96,artifactDataTime*30000/1000);
    artifactData.t = zeros(3000,1);
    artifactDataIndex = 1;
     
    for i=1:numel(fileList)
        %% load file
        disp(['working on:'])
        disp(fileList(i).name)
        
        cdsTemp=commonDataStructure();

        cdsTemp.file2cds([folderpath,fileList(i).name],inputData.ranBy,inputData.array1,inputData.monkey,inputData.lab,inputData.task,inputData.mapFile);
        
        
        spikeWaves = zeros(10000,lengthWave);
        spikeTimes = zeros(10000,1);
        spikeChan = zeros(10000,1);    
        spikeNum = 1;
        
        %% go through channels and get thresholdCrossings
        % write this later when I care?
        
        
        %% get spike data from ainp15
        
        data = -1*cdsTemp.analog{1,1}.ainp15;
        % filter data
        [b,a] = butter(4,[250,5000]/(30000/2),'bandpass');
        data = filtfilt(b,a,data);
        
        % get threshold crossings
        thresholdCrossings = find(data<=-4);

        % remove potential artifacts and too close to beginning/end
        % of stim data
        crossingsMask = ones(numel(thresholdCrossings),1);
        for cross = 1:numel(thresholdCrossings)
            if(~(thresholdCrossings(cross)+postOffset <= numel(data) && thresholdCrossings(cross)-preOffset > 0))
                crossingsMask(cross) = 0;
            end
        end
        thresholdCrossings = thresholdCrossings(crossingsMask(:) == 1);

        % remove chains -- find largest spot
        idx = 2;
        chain = [1];
        crossingsKeep = [];
        while idx <= numel(thresholdCrossings)
            if(thresholdCrossings(idx) == thresholdCrossings(idx-1)+1) % store in chain
                chain = [chain;idx];
            elseif(~isempty(chain)) % broke a chain, store minidx, update idx, empty chain
                [~,minIdx] = min(data(thresholdCrossings(chain)));
                if(isempty(crossingsKeep))
                    crossingsKeep = [thresholdCrossings(minIdx+chain(1)-1)];
                else
                    crossingsKeep = [crossingsKeep;thresholdCrossings(minIdx+chain(1)-1)];
                end
                chain = [idx];
            end
            idx = idx+1;
        end
        if(numel(thresholdCrossings) > 0)
            thresholdCrossings = [crossingsKeep;thresholdCrossings(end)];
        end

        % go through and weed out ones that are too close to each other
        % prioritize backwards in time
        crossingsMask = ones(numel(thresholdCrossings),1);
        for cross = numel(thresholdCrossings):-1:2
            if(crossingsMask(cross) == 1) % check time beforehand to see if one is too close
                crossCheck = cross-1;
                while crossCheck >= 1 && thresholdCrossings(crossCheck) >= thresholdCrossings(cross) - max(preOffset,postOffset)
                    crossingsMask(crossCheck) = 0;
                    crossCheck = crossCheck-1;
                end
            end
        end  
        thresholdCrossings = thresholdCrossings(crossingsMask(:) == 1);

        % store thresholdCrossing data
        for cross = 1:numel(thresholdCrossings)
            spikeTimes(spikeNum) = cdsTemp.analog{1,1}.t(thresholdCrossings(cross)); % this is in seconds

            spikeChan(spikeNum) = 1;
            % check if too close to beginning or end
            if(thresholdCrossings(cross)+postOffset > numel(data))
                numZerosPad = (preOffset+postOffset+1) - numel(data(thresholdCrossings(cross)-preOffset:end));
                spikeWaves(spikeNum,:) = [data(thresholdCrossings(cross)-preOffset:end);zeros(numZerosPad,1)];
            elseif(thresholdCrossings(cross)-preOffset < 0)
                numZerosPad = (preOffset+postOffset+1) - numel(data(1:thresholdCrossings(cross)+postOffset));
                spikeWaves(spikeNum,:) = [zeros(numZerosPad,1);data(thresholdCrossings(cross)-preOffset:end)];
            else
                spikeWaves(spikeNum,:) = data(thresholdCrossings(cross)-preOffset:thresholdCrossings(cross)+postOffset);
            end
            spikeNum = spikeNum + 1;

            if(spikeNum > 0.67*numel(spikeTimes))
                spikeTimes = [spikeTimes;zeros(1000,1)];
                spikeWaves = [spikeWaves; zeros(1000,lengthWave)];
                spikeChan = [spikeChan; zeros(1000,1)];
            end
        end
        
        
    end
    
    
    %% write nev file
    disp('writing nev file')
    spikeTimes = spikeTimes(1:spikeNum-1,1);
    spikeWaves = spikeWaves(1:spikeNum-1,:);
    spikeChan = spikeChan(1:spikeNum-1,:);
    [spikeTimes,sortOrder] = sort(spikeTimes);
    spikeWaves = spikeWaves(sortOrder,:);
    spikeChan = spikeChan(sortOrder,:);

    % remove spikeWaves that are too large
    spikeWaves = spikeWaves/100*1000*4;
    maskWaves = max(abs(spikeWaves),[],2) < 1000;
    spikeWaves = spikeWaves(maskWaves,:);
    spikeTimes = spikeTimes(maskWaves);
    spikeChan = spikeChan(maskWaves);
    
    nevData.ts = spikeTimes;
    nevData.waveforms = spikeWaves(:,:);
    nevData.elec = spikeChan(:,:);

    packetWidth = 104;
    filename = strcat(fileList(1).name(1:end-4),'_spikeInformation');
    mapFilename = inputData.mapFile(8:end);
    comments = '';
    writeNEV(nevData, packetWidth, filename, mapFilename, comments )
    
end