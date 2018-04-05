function [arrayData,report] = labelResponse(arrayData,opts)

    %% setup options
    opts = configureOpts(opts);
    


    %% for each unit, compute mean and SD from baseline. Then check for excitatory
    % and inhibitory
    for chan = 1:size(arrayData{1}.spikeTrialTimes,1)
        for wave = 1:size(arrayData{1}.spikeTrialTimes,2)
            report.excitatoryUnits{chan,wave} = [];
            report.inhibitoryUnits{chan,wave} = [];
            
            for arrayIdx = 1:numel(arrayData)
                % initialize storage
                arrayData{arrayIdx}.isExcitatory{chan,wave} = 0;
                arrayData{arrayIdx}.excitatoryLatency{chan,wave} = [1,1,1]; % onset, peak, end
                arrayData{arrayIdx}.isInhibitory{chan,wave} = 0;
                arrayData{arrayIdx}.inhibitoryLatency{chan,wave} = [1,1]; % onset, end

                % find bin edges
                baseline_binEdgePre = max(find(arrayData{arrayIdx}.bE{chan,wave} <= opts.BASELINE_PRE_TIME*1000));
                baseline_binEdgePost = max(find(arrayData{arrayIdx}.bE{chan,wave} <= opts.BASELINE_POST_TIME*1000)) - 1; % one more bin edge than bin count
                minBin = max(find(arrayData{arrayIdx}.bE{chan,wave} <= opts.MIN_TIME*1000));
                maxExcitatoryBin = max(find(arrayData{arrayIdx}.bE{chan,wave} <= opts.MAX_EXCITATION_TIME*1000));
                maxInhibitoryBinEnd = max(find(arrayData{arrayIdx}.bE{chan,wave} <= opts.MAX_INHIBITION_TIME_END*1000));
                maxInhibitoryBinStart = max(find(arrayData{arrayIdx}.bE{chan,wave} <= opts.MAX_INHIBITION_TIME_START*1000));
               
                % compute mean and SD
                mean_BC = mean(arrayData{arrayIdx}.bC{chan,wave}(baseline_binEdgePre:baseline_binEdgePost));
                SD_BC = std(arrayData{arrayIdx}.bC{chan,wave}(baseline_binEdgePre:baseline_binEdgePost));

                % get windows
                if(opts.BRUTE_FORCE_MAX_WINDOW) % finds the window that results in the maximum response
                    maxExcitatoryResponse = -10000;
                    minInhibitoryResponse = 100000;
                    arrayData{arrayIdx}.excitatoryLatency{chan,wave} = [];
                    % brute force all possible windows
                    for binStart = minBin:maxExcitatoryBin-1
                        for binEnd = binStart+opts.MIN_WINDOW_SIZE:maxExcitatoryBin
                            excitatoryResponse = sum(arrayData{arrayIdx}.bC{chan,wave}(binStart:binEnd)) - (opts.BRUTE_FORCE_EXCITE_GAMMA)*mean_BC*(binEnd-binStart); % consider penalizing for window sie
                            if(excitatoryResponse > maxExcitatoryResponse)
                                maxExcitatoryResponse = excitatoryResponse;
                                [~,peakIdx] = max(arrayData{arrayIdx}.bC{chan,wave}(binStart:binEnd));
                                peakIdx = peakIdx+binStart-1;
                                arrayData{arrayIdx}.excitatoryLatency{chan,wave} = [binStart,peakIdx,binEnd];
                            end
                        end
                    end
                    if(sum(arrayData{arrayIdx}.bC{chan,wave}(peakIdx:peakIdx+opts.NUM_EXCITATORY_BINS-1) > mean_BC+2*SD_BC) == opts.NUM_EXCITATORY_BINS)
                        arrayData{arrayIdx}.isExcitatory{chan,wave} = 1;
                        report.excitatoryUnits{chan,wave}(end+1) = arrayIdx;
                    end
                    
                    for binStart = minBin:maxInhibitoryBinStart-1
                        for binEnd = binStart+opts.MIN_WINDOW_SIZE:maxInhibitoryBinEnd
                            inhibitoryResponse = sum(arrayData{arrayIdx}.bC{chan,wave}(binStart:binEnd)) - opts.BRUTE_FORCE_INHIB_GAMMA*(binEnd-binStart);
                            if(inhibitoryResponse < minInhibitoryResponse)
                                minInhibitoryResponse = inhibitoryResponse;
                                arrayData{arrayIdx}.inhibitoryLatency{chan,wave} = [binStart,binEnd];
                            end
                        end
                    end
                    binStart = arrayData{arrayIdx}.inhibitoryLatency{chan,wave};
                    if(sum(arrayData{arrayIdx}.bC{chan,wave}(binStart:binStart+opts.NUM_INHIBITORY_BINS-1) < mean_BC-SD_BC) == opts.NUM_INHIBITORY_BINS)
                        arrayData{arrayIdx}.isInhibitory{chan,wave} = 1;
                        report.inhibitoryUnits{chan,wave}(end+1) = arrayIdx;
                    end
                else % count consecutive peaks, pick first instance. This is not great, but works in general
                    % three consecutive bins above mean + 2*SD = excitatory
                    flag_excitatory = 0;
                    for binIdx = minBin:maxExcitatoryBin-opts.NUM_EXCITATORY_BINS
                        if(~flag_excitatory && sum(arrayData{arrayIdx}.bC{chan,wave}(binIdx:binIdx+opts.NUM_EXCITATORY_BINS-1) > mean_BC+2*SD_BC) == opts.NUM_EXCITATORY_BINS)
                            arrayData{arrayIdx}.isExcitatory{chan,wave} = 1;
                            arrayData{arrayIdx}.excitatoryLatency{chan,wave}(1) = binIdx;
                            arrayData{arrayIdx}.excitatoryLatency{chan,wave}(2) = binIdx;
                            flag_excitatory = 1;
                            report.excitatoryUnits{chan,wave}(end+1) = arrayIdx;
                        elseif(flag_excitatory && arrayData{arrayIdx}.bC{chan,wave}(binIdx) > arrayData{arrayIdx}.bC{chan,wave}(arrayData{arrayIdx}.excitatoryLatency{chan,wave}(2)))
                            % find peak
                            arrayData{arrayIdx}.excitatoryLatency{chan,wave}(2) = binIdx;
                        elseif(flag_excitatory && sum(arrayData{arrayIdx}.bC{chan,wave}(binIdx:binIdx+opts.NUM_EXCITATORY_BINS-1) < mean_BC+2*SD_BC) < opts.MIN_EXCITATORY_BINS)
                            % find end of excitatory (break out of for loop)
                            arrayData{arrayIdx}.excitatoryLatency{chan,wave}(3) = binIdx;
                            break;
                        end
                    end

                    % five consecutive bins below mean - SD = inhibitory
                    flag_inhibitory = 0;
                    for binIdx = minBin:maxInhibitoryBinEnd-opts.NUM_INHIBITORY_BINS
                        if(~flag_inhibitory && sum(arrayData{arrayIdx}.bC{chan,wave}(binIdx:binIdx+opts.NUM_INHIBITORY_BINS-1) < mean_BC-SD_BC) == opts.NUM_INHIBITORY_BINS)
                            arrayData{arrayIdx}.isInhibitory{chan,wave} = 1;
                            arrayData{arrayIdx}.inhibitoryLatency{chan,wave}(1) = binIdx;
                            arrayData{arrayIdx}.inhibitoryLatency{chan,wave}(2) = binIdx;
                            flag_inhibitory = 1;
                            report.inhibitoryUnits{chan,wave}(end+1) = arrayIdx;
                        elseif(flag_inhibitory && sum(arrayData{arrayIdx}.bC{chan,wave}(binIdx:binIdx+opts.NUM_INHIBITORY_BINS-1) < mean_BC-SD_BC) < opts.MIN_INHIBITORY_BINS)
                            % find end of excitatory (break out of for loop)
                            arrayData{arrayIdx}.inhibitoryLatency{chan,wave}(2) = binIdx;
                            break;
                        end
                    end

                end
                % set binIdxs to times
                arrayData{arrayIdx}.excitatoryLatency{chan,wave} = arrayData{arrayIdx}.bE{chan,wave}(arrayData{arrayIdx}.excitatoryLatency{chan,wave}) + mode(diff(arrayData{arrayIdx}.bE{chan,wave}))/2;
                arrayData{arrayIdx}.inhibitoryLatency{chan,wave} = arrayData{arrayIdx}.bE{chan,wave}(arrayData{arrayIdx}.inhibitoryLatency{chan,wave}) + mode(diff(arrayData{arrayIdx}.bE{chan,wave}))/2;
                
            end
            
            report.excitatoryUnits{chan,wave} = unique(report.excitatoryUnits{chan,wave});
            report.inhibitoryUnits{chan,wave} = unique(report.inhibitoryUnits{chan,wave});
        end
    end


    
end


function [opts] = configureOpts(optsInput)

    opts.BASELINE_PRE_TIME = -19/1000;
    opts.BASELINE_POST_TIME = -2/1000;

    opts.MIN_TIME = 1.2/1000;

    opts.MAX_INHIBITION_TIME_END = 75/1000;
    opts.MAX_INHIBITION_TIME_START = 20/1000;
    opts.MAX_EXCITATION_TIME = 20/1000;

    opts.NUM_EXCITATORY_BINS = 3;
    opts.NUM_INHIBITORY_BINS = 5;
    opts.MIN_INHIBITORY_BINS = 3;
    opts.MIN_EXCITATORY_BINS = 3;
    
    opts.MIN_WINDOW_SIZE = 3;
    opts.BRUTE_FORCE_MAX_WINDOW = 0;
    opts.BRUTE_FORCE_EXCITE_GAMMA = 3;
    opts.BRUTE_FORCE_INHIB_GAMMA = 3;
    
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