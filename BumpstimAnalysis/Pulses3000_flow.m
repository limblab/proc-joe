folderpath = 'C:\Users\Joseph\Desktop\Lab\Data\StimArtifact\Han_20170621_3000stimuli\';
pwd = cd;
cd(folderpath)
folderList = dir('chan*');
for fold = 1:numel(folderList)
    disp(fold)
    filename = folderList(fold).name;
    nevFilename = strcat('Han_20170516_',filename,'500HzFilter');
    numCrossings = 1;
    [b,a] = butter(6,500/(30000/2),'high');
    for garb = 1:3
        disp(garb)
        load(strcat(folderpath,folderList(fold).name,filesep,'artifactData.mat'));
        load(strcat(folderpath,folderList(fold).name,filesep,'chList.mat'));
        load(strcat(folderpath,folderList(fold).name,filesep,'eList.mat'));
        load(strcat(folderpath,folderList(fold).name,filesep,'posList.mat'));
        outputData.artifactData = artifactData; clear artifactData;
        outputData.chList = chList; clear chList;
        outputData.eList = eList; clear eList;
        outputData.posList = posList; clear posList;
        %
        if(garb==1)
            outputData.artifactData.artifact = outputData.artifactData.artifact(:,1:1000,:);
            outputData.artifactData.stimOn = outputData.artifactData.stimOn(1:1000);
        elseif(garb==2)
            outputData.artifactData.artifact = outputData.artifactData.artifact(:,1001:2000,:);
            outputData.artifactData.stimOn = outputData.artifactData.stimOn(1001:2000);
        else
            outputData.artifactData.artifact = outputData.artifactData.artifact(:,2001:end,:);
            outputData.artifactData.stimOn = outputData.artifactData.stimOn(2001:end);
        end

        % filter
        for ch = 1:96
            outputData.artifactData.artifact(ch,:,:) = squeeze(fliplr(filter(b,a,fliplr(squeeze(outputData.artifactData.artifact(ch,:,:)))')'));
        end

        % get threshold
        thresholdAll = zeros(96,1);
        idx10 = 30*30;
        numSamples = (size(outputData.artifactData.artifact,3)-idx10+1)*size(outputData.artifactData.artifact,2);
        for ch = 1:96
            for stimuli = 1:size(outputData.artifactData.artifact,2)
                thresholdAll(ch,1) = thresholdAll(ch,1) + sum(outputData.artifactData.artifact(ch,stimuli,idx10:end).^2)/numSamples;
            end
            thresholdAll(ch,1) = 3*sqrt(thresholdAll(ch,1));
        end

        % threshold and save waveforms (nevData)
        idx5 = 5*30;
        idx6 = 6*30;
        idx1 = 1*30;
        preOffset = 27;
        postOffset = 20;
        lengthWave = postOffset + preOffset + 1; % plus one for the 0 idx
        % ts = zeros(1000,1);
        % waves = zeros(1000,lengthWave);
        % chan = zeros(1000,1);
        % numCrossings = 1;

        for stim = 1:size(outputData.artifactData.artifact,2)
            for elec = 1:size(outputData.artifactData.artifact,1)
                stimData = squeeze(outputData.artifactData.artifact(elec,stim,idx5:end));
                stimData(1:idx1) = 0; % blank the artifact data
                threshold = abs(thresholdAll(elec));
                thresholdCrossings = find(stimData(:)>threshold);

                % remove chains -- find min in chain and keep that        
                idx = 2;
                chain = [1];
                crossingsKeep = [];
                while idx <= numel(thresholdCrossings)
                    if(thresholdCrossings(idx) == thresholdCrossings(idx-1)+1) % store in chain
                        chain = [chain;idx];
                    elseif(~isempty(chain)) % broke a chain, store minidx, update idx, empty chain
                        [~,maxIdx] = max(stimData(thresholdCrossings(chain)));
                        if(isempty(crossingsKeep))
                            crossingsKeep = [thresholdCrossings(maxIdx+chain(1)-1)];
                        else
                            crossingsKeep = [crossingsKeep;thresholdCrossings(maxIdx+chain(1)-1)];
                        end
                        chain = [idx];
                    end
                    idx = idx+1;
                end
                if(numel(thresholdCrossings) > 0)
                    thresholdCrossings = [crossingsKeep;thresholdCrossings(end)];
                end
                % remove potential artifacts and too close to beginning
                crossingsMask = ones(numel(thresholdCrossings),1);
                for cross = 1:numel(thresholdCrossings)
                    if(stimData(thresholdCrossings(cross)) > 1000 || ...
                            ~(thresholdCrossings(cross)+postOffset <= numel(stimData) && thresholdCrossings(cross)-preOffset > 0))
                        crossingsMask(cross) = 0;
                    end
                end
                thresholdCrossings = thresholdCrossings(crossingsMask(:) == 1);

                % go through and weed out ones that are too close - prioritize
                % backwards in time
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
                    ts(numCrossings) = outputData.artifactData.stimOn(stim)/30000 + thresholdCrossings(cross)/30000 - 1/30000; % this is in seconds
                    chan(numCrossings) = elec;
                    waves(numCrossings,:) = stimData(thresholdCrossings(cross)-preOffset:thresholdCrossings(cross)+postOffset);
                    numCrossings = numCrossings + 1;

                    if(numCrossings > 0.67*numel(ts))
                        ts = [ts;zeros(1000,1)];
                        waves = [waves; zeros(1000,lengthWave)];
                        chan = [chan; zeros(1000,1)];
                    end
                end
            end 
        end


    end % garb

    % save nevData
    nevData.ts = ts(1:numCrossings);
    nevData.waveforms = waves(1:numCrossings,:);
    nevData.elec = chan(1:numCrossings);
    packetWidth = (postOffset+preOffset+1)*2 + 8; % num bytes in data = length of data*2 + 8 bytes
    mapFilename = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    comments = '';
    [nevData.ts,sortTs] = sort(nevData.ts);
    nevData.waveforms = nevData.waveforms(sortTs,:);
    nevData.elec = nevData.elec(sortTs,:);
    save(nevFilename,'nevData');
    disp('done saving mat files')
end

%% merge .mat files
fileList  = dir('*Filter.mat');
for f = 1:numel(fileList)
    load(fileList(f).name);
    
    if(f==1)
        nevDataAll.ts = nevData.ts;
        nevDataAll.waveforms = nevData.waveforms;
        nevDataAll.elec = nevData.elec;
    else
        nevDataAll.ts(end+1:end+numel(nevData.ts)) = nevData.ts + (f-1)*1000;
        nevDataAll.waveforms(end+1:end+numel(nevData.ts),:) = nevData.waveforms;
        nevDataAll.elec(end+1:end+numel(nevData.ts)) = nevData.elec;
    end
end

%%
packetWidth = (postOffset+preOffset+1)*2 + 8; % num bytes in data = length of data*2 + 8 bytes
mapFilename = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
comments = '';
nevFilename = 'Han_20170516_all_500HzFilter-positiveThreshold';
disp('writing NEV file')
writeNEV(nevDataAll, packetWidth, nevFilename, mapFilename, comments);
disp('done writing NEV file')
