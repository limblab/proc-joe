function [ common ] = getAllNeuronRates(filepath,files,neurons,varargin )
% this function gets all of the firing rate histograms for the neurons in
% neurons for all files in files.
binSize = 0.05;
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'binSize'
            binSize = varargin{i+1};
    end
end

% get electrode/channel for each file/neuron in table of neurons
[common.muscles, common.channels, common.IDs, common.electrodes, common.neuronNumbers] = ... 
    getElectrodeChannels(neurons,files, filepath);

% remove neurons that are in common multiple times
common = removeDuplicateNeurons(common);

% for each file
common.muscleNames = {};
for i = 1:length(files)
    % get the muscle name for the file -- this is for data storage purposes
    idxUnderScore = strfind(files(i).name,'_');
    muscleName = files(i).name(idxUnderScore(3)+1:idxUnderScore(4)-1);

    % load cds from filename
    load(strcat(filepath,files(i).name));

    % for every neuron in common, find neuron number correspoding to
    % channel and ID and elec, then generatePESTH(cds,nn,'noPlots',1). This
    % will output the counts/edges for that neuron. Store this as
    % common.<muscle name>.<counts or edges>. 
    common.(genvarname(muscleName)).counts = {};
    common.(genvarname(muscleName)).edges = {};
    common.(genvarname(muscleName)).bootstrapMean = [];
    common.(genvarname(muscleName)).bootstrapBounds = [];
    common.(genvarname(muscleName)).meanCount = [];
    common.(genvarname(muscleName)).stdCount = [];
    common.muscleNames{end+1} = muscleName;

    % find sequenceTimes and eventTimes once
    GTOstim = 0;
    if(size(cds.analog,1) == 0)
        GTOstim = 1;
    end
    if(~exist('sequenceTimes') || ~exist('eventTimes') || ~exist('stimState'))
        [stimState,] = determineStimTiming(cds, GTOstim, 0);
        [sequenceTimes, eventTimes] = getSequenceTimes(cds, stimState,GTOstim,0);
    end
    common.(genvarname(muscleName)).eventTimes = eventTimes;
    common.(genvarname(muscleName)).sequenceTimes = sequenceTimes;
    common.(genvarname(muscleName)).stimState = stimState;
    
    for j = 1:length(common.muscles)
        nn = findNeuronNumber(cds, common.channels(j), common.IDs(j), common.electrodes{j});
        if(nn == -1) % generate PSTH to get correct sizes, fill in counts with zeros
            [bE, bC] = generatePESTH(cds,1,'useRate',1,'noPlots',1,'sequenceTimes',sequenceTimes,'eventTimes',eventTimes,...
                'stimState',stimState,'binSize',binSize);
            common.(genvarname(muscleName)).counts{end+1,1} = zeros(size(bC));
            common.(genvarname(muscleName)).edges{end+1,1} = bE;
            common.(genvarname(muscleName)).bootstrapMean(end+1,1) = 0;
            common.(genvarname(muscleName)).bootstrapBounds(end+1,:) = [0,0];
            common.(genvarname(muscleName)).meanCount(end+1,1) = 0;
            common.(genvarname(muscleName)).stdCount(end+1,1) = 0;
        else
            [bE, bC] = generatePESTH(cds,nn,'useRate',1,'noPlots',1,'sequenceTimes',sequenceTimes,'eventTimes',eventTimes,...
                'stimState',stimState,'binSize',binSize);
            
            common.(genvarname(muscleName)).counts{end+1,1} = bC;
            common.(genvarname(muscleName)).edges{end+1,1} = bE;
            [common.(genvarname(muscleName)).bootstrapMean(end+1,1), common.(genvarname(muscleName)).bootstrapBounds(end+1,:)] = ...
                bootstrapConfidenceInterval(cds,nn,eventTimes(1),binSize);
            common.(genvarname(muscleName)).meanCount(end+1,1) = mean(bC);
            common.(genvarname(muscleName)).stdCount(end+1,1) = std(bC);
        end
    end
end

end

