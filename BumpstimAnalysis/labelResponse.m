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
                
                maxExcitatoryBinEnd = max(find(arrayData{arrayIdx}.bE{chan,wave} <= opts.MAX_EXCITATION_TIME_END*1000));
                maxExcitatoryBinStart = max(find(arrayData{arrayIdx}.bE{chan,wave} <= opts.MAX_EXCITATION_TIME_START*1000));
                maxInhibitoryBinEnd = max(find(arrayData{arrayIdx}.bE{chan,wave} <= opts.MAX_INHIBITION_TIME_END*1000));
                maxInhibitoryBinStart = max(find(arrayData{arrayIdx}.bE{chan,wave} <= opts.MAX_INHIBITION_TIME_START*1000));
               
                % compute mean and SD
                mean_BC = mean(arrayData{arrayIdx}.bC{chan,wave}(baseline_binEdgePre:baseline_binEdgePost));
                SD_BC = std(arrayData{arrayIdx}.bC{chan,wave}(baseline_binEdgePre:baseline_binEdgePost));

                %% find peak bin, then expand 
                [~,largestBin] = max(arrayData{arrayIdx}.bC{chan,wave}(minBin:maxExcitatoryBinStart));
                largestBin = largestBin + minBin - 1; % shift the bin idx
                arrayData{arrayIdx}.excitatoryLatency{chan,wave} = [largestBin,largestBin,largestBin];
                
                % expand window backwards 
                for idx = largestBin:-1:minBin
                    if(arrayData{arrayIdx}.bC{chan,wave}(idx) > mean_BC + SD_BC)
                        arrayData{arrayIdx}.excitatoryLatency{chan,wave}(1) = idx;
                    else
                        break;
                    end
                end
                    
                % expand window forwards
                for idx = largestBin:1:maxExcitatoryBinEnd
                    if(arrayData{arrayIdx}.bC{chan,wave}(idx) > mean_BC + SD_BC)
                        arrayData{arrayIdx}.excitatoryLatency{chan,wave}(3) = idx;
                    else
                        break;
                    end
                end
                % label based on peak
                if(sum(arrayData{arrayIdx}.bC{chan,wave}(largestBin-1:largestBin+1) > mean_BC + SD_BC) >= 3)
                    arrayData{arrayIdx}.isExcitatory{chan,wave} = 1;
                    report.excitatoryUnits{chan,wave}(end+1) = arrayIdx;
                end
                
                %% find smallest bin, then expand
                [~,smallestBin] = min(arrayData{arrayIdx}.bC{chan,wave}(maxInhibitoryBinStart:maxInhibitoryBinEnd));
                smallestBin = smallestBin + maxInhibitoryBinStart - 1; % shift the bin idx
                arrayData{arrayIdx}.inhibitoryLatency{chan,wave} = [smallestBin,smallestBin];
                
                % expand window backwards 
                for idx = smallestBin:-1:minBin
                    if(arrayData{arrayIdx}.bC{chan,wave}(idx) < mean_BC - 0.5*SD_BC)
                        arrayData{arrayIdx}.inhibitoryLatency{chan,wave}(1) = idx;
                    else
                        break;
                    end
                end
                    
                % expand window forwards
                for idx = smallestBin:1:maxInhibitoryBinEnd
                    if(arrayData{arrayIdx}.bC{chan,wave}(idx) < mean_BC - 0.5*SD_BC)
                        arrayData{arrayIdx}.inhibitoryLatency{chan,wave}(2) = idx;
                    else
                        break;
                    end
                end
                % label based on trough
                if(sum(arrayData{arrayIdx}.bC{chan,wave}(smallestBin-2:smallestBin+2) < mean_BC - SD_BC) >= 4)
                    arrayData{arrayIdx}.isInhibitory{chan,wave} = 1;
                    report.inhibitoryUnits{chan,wave}(end+1) = arrayIdx;
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
    opts.MAX_TIME = 75/1000;
    
    opts.MAX_INHIBITION_TIME_END = 75/1000;
    opts.MAX_INHIBITION_TIME_START = 15/1000;
    opts.MAX_EXCITATION_TIME_START = 5/1000;
    opts.MAX_EXCITATION_TIME_END = 7/1000;
    
  
    
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