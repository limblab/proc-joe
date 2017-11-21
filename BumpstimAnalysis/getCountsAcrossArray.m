function [outputData] = getCountsAcrossArray(cds,opts)

    % configure opts
    opts = configureOpts(opts);
    
    % check for bad inputs
    if(strcat(opts.MAP_FILENAME,'')==1)
        error('no map file');
    end
    
    try
        mapData=loadMapFile(opts.MAP_FILENAME);
    catch
        error('bad map file')
    end
    
    % establish outputData
    outputData = [];
    
    % loop through each neuron in cds and get bC and bE for each waveform
    % type and channel requested
    posList = [mapData.row,mapData.col];
    outputDataIdx = 1;
    for nn = 1:size(cds.units,2)
        if(cds.units(nn).ID ~= 0 && cds.units(nn).ID ~= 255) % ignore unsorted and invalid
            posIdx=find(mapData.chan==cds.units(nn).chan);
            eRow=posList(posIdx,1);
            eRow = 11 - eRow;
            eCol=posList(posIdx,2);
            
            outputData{outputDataIdx}.neuronNumber = nn;
            outputData{outputDataIdx}.row = eRow;
            outputData{outputDataIdx}.col = eCol;
            
            if(opts.BIN_DATA == 1) % will get both spike and bin data, may not store spike data though
                % get bin counts and edges using plotPSTHStim
                [data] = plotPSTHStim(cds,nn,'binSize',opts.BIN_SIZE,'waveformTypes',opts.WAVEFORM_IDX,...
                    'chans',opts.CHANNEL_IDX,'preTime',opts.PRE_TIME,'postTime',opts.POST_TIME,'noPlot',1);

                % put counts and row/col info into outputData
                outputData{outputDataIdx}.parameters = cds.waveforms.parameters;
                outputData{outputDataIdx}.chanSent = unique(cds.waveforms.chanSent);
                
                outputData{outputDataIdx}.bC = data.bC;
                outputData{outputDataIdx}.bE = data.bE;
                if(opts.SPIKE_DATA)
                    outputData{outputDataIdx}.spikeTrialTimes = data.spikeTrialTimes;
                    outputData{outputDataIdx}.spikeTrueTimes = data.spikeTrueTimes;
                    outputData{outputDataIdx}.stimData = data.stimData;
                    outputData{outputDataIdx}.numStims = data.numStims;
                end

            elseif(opts.SPIKE_DATA == 1)
                % get bin counts and edges using plotPSTHStim
                [data] = plotRasterStim(cds,nn,'binSize',opts.BIN_SIZE,'waveformTypes',opts.WAVEFORM_IDX,...
                    'chans',opts.CHANNEL_IDX,'preTime',opts.PRE_TIME,'postTime',opts.POST_TIME,'noPlot',1);

                outputData.spikeTimes{outputDataIdx} = data.spikeTimes;
                outputData.stimData{outputDataIdx} = data.stimData;
 
            end
            outputDataIdx = outputDataIdx + 1;
        end
    end
    
end

function [opts] = configureOpts(optsInput)

    opts.MAP_FILENAME = '';
    opts.PRE_TIME = 10/1000;
    opts.POST_TIME = 30/1000;
    opts.BIN_SIZE = 0.0002;
    opts.WAVEFORM_IDX = 1;
    opts.CHANNEL_IDX = 1;
    opts.BIN_DATA = 1;
    opts.SPIKE_DATA = 0;
    
    %% check if in opts and optsInput, overwrite if so
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