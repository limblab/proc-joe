function [  ] = compareAcrossFiles( common, neuronNumberCommon, filepaths,fileprefixes )

for i = 1:numel(filepaths) % for each file to compare to
    % get all files that matter
    filepath = filepaths{i}; fileprefix = fileprefixes{i};
    files = dir(strcat(filepath, fileprefix, '*_cds.mat'));
    
    % remove files with muscles that aren't in common and significant
    j = 1;
    while j <= length(files)
        % get the muscle name for the file -- this is for data storage purposes
        idxUnderScore = strfind(files(j).name,'_');
        muscleName = files(j).name(idxUnderScore(3)+1:idxUnderScore(4)-1);
        muscleIdx=(strfind(common.neuronLabels{neuronNumberCommon},muscleName));
        keep = 0;
        for k = 1:length(muscleIdx)
            if(~isempty(muscleIdx{k}))
                keep = 1;
            end
        end
        if(keep == 0)
            files(j) = [];
            j = j-1;
        end
        j = j+1;
    end
    disp(numel(files))
    % go through each file in files
    for j = 1:length(files)
        % load cds from that file
        load(strcat(filepath,files(j).name));
        
        % find neuron number in the cds
        neuronNumbersCDS = [];
        flagStop = 0;
        idx = 0;
        while(flagStop == 0)
            nn = findNeuronNumber(cds,common.channels(neuronNumberCommon),idx,common.electrodes{neuronNumberCommon});
            if(nn == -1) % not found
                flagStop = 1;
            else
                neuronNumbersCDS(end+1,1) = nn;
            end
            idx = idx+1;
        end
        
        % plot PSTH for neuron in common and neuronNumbersCDS
        figure(); 
        [p,] = numSubplots(1+numel(neuronNumbersCDS));
        % plot common neuron
%         subplot(p(1),p(2),1)
        idxUnderScore = strfind(files(j).name,'_');
        muscleName = files(j).name(idxUnderScore(3)+1:idxUnderScore(4)-1);
        plotPSTHCommon(common,neuronNumberCommon,muscleName,'confidenceInterval',1,'useRate',1,'figMake',0)
        title(common.fileprefix,'Interpreter','none');
        
        % plot other neurons
        binEdges = common.(muscleName).edges{1};
        binSize = binEdges(2)-binEdges(1);
        for k = 1:length(neuronNumbersCDS)
            figure;
%             subplot(p(1),p(2),k+1)
            generatePESTH(cds,neuronNumbersCDS(k),'useRate',1, 'averageSpikeWaveform',0,...
                'sequenceTimes',sequenceTimes,'eventTimes',eventTimes,'confidenceInterval',1,...
                'binSize',binSize,'binsAbove',0,'figMake',0,'legendMake',0,'highlightBin',1);
            title(fileprefixes(i),'Interpreter','none');
        end
%         suptitle(muscleName)  
        
        % plot average waveforms :( -this is super poorly written btw
        cds2 = cds; % so i can load the cds for the common neuron -- lol
        commonFiles = dir(strcat(common.filepath, common.fileprefix, '*_cds.mat'));
        fileIdx = 0;
        for idx = 1:length(commonFiles)
            idxUnderScore = strfind(commonFiles(idx).name,'_');
            fileMuscleName = commonFiles(idx).name(idxUnderScore(3)+1:idxUnderScore(4)-1);
            if(strcmp(muscleName,fileMuscleName)==1)
                fileIdx = idx;
            end
        end
        load(strcat(common.filepath,commonFiles(fileIdx).name));
        figure();
%         [p,] = numSubplots(1+numel(neuronNumbersCDS));
        % plot average waveform for common neuron
        nn = common.neuronNumbers(neuronNumberCommon);
        subplot(p(1),p(2),1)
        plotAllWaveforms(cds,nn);
        title(common.fileprefix,'Interpreter','none');
        
        % plot average waveform for other neurons
        for k = 1:length(neuronNumbersCDS)
            figure
%             subplot(p(1),p(2),k+1)
            plotAllWaveforms(cds2,neuronNumbersCDS(k));
            title(fileprefixes(i),'Interpreter','none');
        end
%         suptitle(muscleName)  
        
        
        % plot average waveforms 
        figure();
        [p,] = numSubplots(1+numel(neuronNumbersCDS));
        % plot average waveform for common neuron
        subplot(p(1),p(2),1)
        plotAverageWaveform(cds,nn);
        title(common.fileprefix,'Interpreter','none');
        
        % plot average waveform for other neurons
        for k = 1:length(neuronNumbersCDS)
            subplot(p(1),p(2),k+1)
            plotAverageWaveform(cds2,neuronNumbersCDS(k));
            title(fileprefixes(i),'Interpreter','none');
        end
%         suptitle(muscleName)  
        
        % remove cds (not really helpful)
        clear cds;
        clear cds2;
    end
    
end

end

