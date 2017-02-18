function [ muscles, channels, IDs, electrodes, neuronNumbers ] = getElectrodeChannels( neurons, fileNames, filePath)
% this function returns the electrodes and channels for all neurons in the
% neuron table provided, which correspond to a file in fileNames
electrodes = {};
channels = [];
IDs = [];
neuronNumbers = [];
muscles = {};
tableNames = neurons.Properties.VariableNames;

for i = 1:length(fileNames) % for each file
    % load cds
    load(strcat(filePath,fileNames(i).name));
    
    % get the relevant muscle from file
    idxUnderScore = strfind(fileNames(i).name,'_');
    muscleName = fileNames(i).name(idxUnderScore(3)+1:idxUnderScore(4)-1);
    
    % find muscleName in neurons table
    % might be multiple instances with other info in name (e.g. amplitude)
    tableIdx = find(strcmp(tableNames,muscleName));
    
    % for each entry in the table related to the muscle in the file
    for t = tableIdx
        nn = neurons{1,t}; % neuron numbers
        for n = nn
            muscles{end+1,1} = muscleName;
            electrodes{end+1,1} = cds.units(n).label;
            channels(end+1,1) = cds.units(n).chan;
            IDs(end+1,1) = cds.units(n).ID;
            neuronNumbers(end+1,1) = n;
        end
    end
    
    clear cds;
end


end

