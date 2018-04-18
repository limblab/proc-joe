function [outputData] = extractDataAroundStimulationsWholeArray(cds,stimInfo,mapFileName,opts)


    %% configure opts and set default values
    opts = configureOpts(opts);
    outputData = [];
    
    %% load in mapfile
    MAP_DATA = loadMapFile(mapFileName);
    POS_LIST = [MAP_DATA.row,MAP_DATA.col];
    
    %% iterate through units in cds and build a cell array with responses
    for unit = 1:size(cds.units,2)
        if(cds.units(unit).ID ~= 0 && cds.units(unit).ID ~= 255)
            % extract and store data
            opts.NEURON_NUMBER = unit;
            unitData = extractDataAroundStimulations(cds,stimInfo,opts);
            
            outputData{end+1,1} = unitData;
            outputData{end}.CHAN = cds.units(unit).chan;
            outputData{end}.ID = cds.units(unit).ID;
            outputData{end}.NN = unit;
            
            posIdx = find(MAP_DATA.chan==cds.units(unit).chan);
            outputData{end}.ROW = 11 - POS_LIST(posIdx,1);
            outputData{end}.COL = POS_LIST(posIdx,2);
            
        end
    end
    
end


function [opts] = configureOpts(optsInput)

    opts.ALIGN_WAVES = 1;
    opts.STIMULI_RESPONSE = 'all'; %'all', 'responsive' or 'nonresponsive'

    opts.STIM_ELECTRODE = 1;
    opts.STIMULATIONS_PER_TRAIN = 1;
    opts.INITIAL_ARRAY_SIZE = 6000;
    opts.ADDITIONAL_ARRAY_SIZE = ceil(opts.INITIAL_ARRAY_SIZE*1/3);
    opts.STIMULATION_BATCH_SIZE = 2000;
    
    opts.PRE_TIME = -5/1000;
    opts.POST_TIME = 30/1000;
    opts.BIN_SIZE = 0.2/1000;
    
    opts.TIME_AFTER_STIMULATION_WAVEFORMS = 10/1000;
    
    opts.NEURON_NUMBER = 1;
    opts.GET_KIN = 1;
    %% check if in optsSave and optsSaveInput, overwrite if so
    try
        inputFieldnames = fieldnames(optsInput);
        for fn = 1:numel(inputFieldnames)
           if(isfield(opts,inputFieldnames{fn}))
               opts.(inputFieldnames{fn}) = optsInput.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
end