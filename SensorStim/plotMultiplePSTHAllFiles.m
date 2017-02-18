% this script plots all of the neurons in '<monkey>_<date>_neurons.mat' as a 
% PSTH in a temperature plot format across all files (muscles typically) in
% the data folder specified. Neurons are organized by the muscle for which
% they were identified

warning('off');
% file prefix :D
filePath = 'D:\Lab\Data\SensorStim\Han_20170209\';
filePrefix = 'Han_20170209';

% check to see if this work has been done already
if(exist(strcat(filePath,filePrefix,'_allNeuronsAllPSTH.mat'))==0)
    % load in neurons <- table with all neuron numbers for each file
    load(strcat(filePath, filePrefix,'_neurons.mat'));

    % get all file names <- cds names
    fileNames = dir(strcat(filePath, filePrefix, '*_cds.mat'));

    % get electrode/channel for each file/neuron in table of neurons
    [common.muscles, common.channels, common.IDs, common.electrodes, common.neuronNumbers] = ... 
        getElectrodeChannels(neurons,fileNames, filePath);

    % remove neurons that are in common multiple times
    common = removeDuplicateNeurons(common);
    
    % get PSTH bin counts and edges for each neuron in each file
    % use generatePESTH(...,'noPlots',1) to get binCounts and binEdges

    % for each file
    common.muscleNames = {};
    for i = 1:length(fileNames)
        % get the muscle name for the file -- this is for data storage purposes
        idxUnderScore = strfind(fileNames(i).name,'_');
        muscleName = fileNames(i).name(idxUnderScore(3)+1:idxUnderScore(4)-1);

        % load cds from filename
        load(strcat(filePath,fileNames(i).name));

        % for every neuron in common, find neuron number correspoding to
        % channel and ID and elec, then generatePESTH(cds,nn,'noPlots',1). This
        % will output the counts/edges for that neuron. Store this as
        % common.<muscle name>.<counts or edges>. 
        common.(genvarname(muscleName)).counts = {};
        common.(genvarname(muscleName)).edges = {};
        common.muscleNames{end+1} = muscleName;

        % find sequenceTimes and eventTimes once
        GTOstim = 0;
        if(size(cds.analog,1) == 0)
            GTOstim = 1;
        end
        [stimState,] = determineStimTiming(cds, GTOstim, 0);
        [sequenceTimes, eventTimes] = getSequenceTimes(cds, stimState,GTOstim,0);  

        for j = 1:length(common.muscles)
            nn = findNeuronNumber(cds, common.channels(j), common.IDs(j), common.electrodes{j});
            [bE, bC] = generatePESTH(cds,nn,'useRate',1,'noPlots',1,'sequenceTimes',sequenceTimes,'eventTimes',eventTimes);
            common.(genvarname(muscleName)).counts{end+1,1} = bC;
            common.(genvarname(muscleName)).edges{end+1,1} = bE;
        end

    end
    save(strcat(filePath,filePrefix,'_allNeuronsAllPSTH'),'common');
else
    load(strcat(filePath,filePrefix,'_allNeuronsAllPSTH'));
end

% so we have all of the bin counts and edges. We need to make an array of
% spike counts for all neurons.matrix will have #neurons x #totalNumBins
% size. Then we can average and stuff, before plotting.

% make matrix of zeros with correct size after getting correct sizes
numNeurons = length(common.muscles);
totalNumBins = getTotalNumBins(common);
common.neuronCounts = zeros(numNeurons, totalNumBins);

% populate common.neuronCounts with the neuronCounts
common = getAllNeuronCounts(common);

% normalize counts relative to average firing rate/counts for each neuron
meanFiring = mean(common.neuronCounts,2);
common.changeNeuronCounts = common.neuronCounts-repmat(meanFiring,[1, size(common.neuronCounts,2)]);
% common.normalizedNeuronCounts(:,:) = min(common.normalizedNeuronCounts(:,:),4);
common.changeNeuronCounts(:,:) = max(common.changeNeuronCounts(:,:),1);

save(strcat(filePath,filePrefix,'_allNeuronsAllPSTH'),'common');


% % plots a bunch of heat map thingies, not super useful
% for i = 1:length(common.muscles)
%     figure();
%     sI = common.muscleIdx{i,2};
%     eI = common.muscleIdx{i,4};
%     imagesc(common.changeNeuronCounts(:,sI:eI));
%     colorbar
%     colormap(hot)
%     hold on
%     zeroIdx = common.muscleIdx{i,3} - common.muscleIdx{i,2}+1;
%     plot([zeroIdx, zeroIdx],[0,40])
%     pause;
%     close all
% end