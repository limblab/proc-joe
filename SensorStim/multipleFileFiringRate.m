function [ firingRate, files ] = multipleFileFiringRate( filepath, fileprefix, neuronNumber )
% given a neuron number and filepath/prefix, find the firingRate of that
% neuron in all files. This assumes that neurons across files will have the
% same CH # and ID $
warning('off')
%% get list of all filenames with prefix
files = dir([filepath fileprefix '*.mat']);
% if there is a -spikes-s file, remove it from list
i = 1;
while i <= length(files)
    if(~isempty(strfind(files(i).name,'-s')))
        files(i) = []; % removes row i from struct files
        i=i-1;
    elseif(~isempty(strfind(files(i).name,'_cds')))
        files(i) = [];
        i=i-1;
    end
    i=i+1;
end

% firingRate holds the firing rate during stim and not during stim for the given
% neuron in each file
firingRate = zeros(length(files),2);

%% For each filename, get firing rate of neuronNumber provided. 

for i = 1:length(files)
    filename = [filepath files(i).name(1:end-4)];
    GTOstim = ~isempty(strfind(lower(filename),'gtostim'));
    timeAfterGTOStim = 0.5; % in seconds
    
    % load cds or give message that cds does not exist and do not move
    % further
    if(exist([filename '_cds.mat'],'file') > 0)
        load([filename '_cds.mat']);
    else
        disp(['File ' files(i).name ' does not have a cds, skipping file']);
        continue; % starts from beginning of loop with next value of i
    end
    
    % determine timing of stimulation
    [stimState,stimDuration,noStimDuration] = determineStimTiming(cds, GTOstim, timeAfterGTOStim);
    firingFrequency = computeFiringFrequency(cds, stimState, stimDuration, noStimDuration, GTOstim);
    firingRate(i,1) = firingFrequency(neuronNumber,1);
    firingRate(i,2) = firingFrequency(neuronNumber,2);
end

warning('on')

end

