function [ outputData ] = thresholdData(artifact,opts)
%
% perform thresholding to get threshold crossings
% threshold = thresholdMult*rms(data(inputData.presample + inputData.windowSize/2:end));
% find all indices with a threshold crossing. Do this by sliding a window backwards,
% filtering with a high pass at like 1 Hz, and then looking at threshold
% crossings based on rms in that window.

    opts = configureOptions(opts);

    % we need to populate outputData.ts (nx1), outputData.waveforms (nx48), and
    % outputData.elec (nx1)
    outputData.ts = [];
    outputData.waveforms = [];
    outputData.elec = [];
    
    for ch = 1:size(artifact,2)
        % filter artifact data
        
        padArtifact = [squeeze(artifact(:,ch,:)),repmat(squeeze(artifact(:,ch,end)),1,opts.NUM_PAD)];
        filterArtifact = fliplr(filter(opts.B,opts.A,fliplr(padArtifact)')');
        filterArtifact = filterArtifact(:,1:size(filterArtifact,2) - opts.NUM_PAD);
        
        threshold = abs(opts.THRESHOLD_MULT)*rms(filterArtifact(1,end-opts.THRESHOLD_WINDOW_LENGTH+1:end));
        outputData.threshold(ch) = threshold;
        
        % get crossings on each artifact
        for art = 1:size(filterArtifact,1)
            crossingIdx = find(squeeze(filterArtifact(art,:)) > threshold);
            
            % process crossings
            
            % remove too close to beginning/end of file and too large
            crossingsMask = ones(numel(crossingIdx),1);
            for cross = 1:numel(crossingIdx)
                if(filterArtifact(art,crossingIdx(cross)) > opts.MAX_AMPLITUDE || ...
                        ~(crossingIdx(cross)+opts.POST_OFFSET <= size(filterArtifact,2) && crossingIdx(cross)-opts.PRE_OFFSET > 0))
                    crossingsMask(cross) = 0;
                end
            end
            crossingIdx = crossingIdx(crossingsMask(:) == 1);
            
            % remove indexes with long strings of crossings, keep max
            idx = 2;
            chain = [1];
            crossingsKeep = [];
            while idx < numel(crossingIdx)
                if(crossingIdx(idx) == crossingIdx(idx-1) + 1)
                    chain = [chain;idx];
                elseif(~isempty(chain))
                    [~,maxIdx] = max(squeeze(filterArtifact(art,crossingIdx(chain))));
                    if(isempty(crossingsKeep))
                        crossingsKeep = [crossingIdx(maxIdx+chain(1)-1)];
                    else
                        crossingsKeep = [crossingsKeep;crossingIdx(maxIdx+chain(1)-1)];
                    end
                    chain = [idx];
                end
                idx = idx + 1;
            end
            if(numel(crossingIdx) > 0)
                crossingIdx = [crossingsKeep;crossingIdx(end)];
            end
            % weed out ones that are too close to each other
            crossingsMask = ones(numel(crossingIdx),1);
            for cross = numel(crossingIdx):-1:2
                if(crossingsMask(cross) == 1) % check time beforehand to see if one is too close
                    crossCheck = cross-1;
                    while crossCheck >= 1 && crossingIdx(crossCheck) >= crossingIdx(cross) - max(opts.PRE_OFFSET,opts.POST_OFFSET)
                        crossingsMask(crossCheck) = 0;
                        crossCheck = crossCheck-1;
                    end
                end
            end  
            crossingIdx = crossingIdx(crossingsMask(:) == 1);

            % store data
            for cross = 1:numel(crossingIdx)
                outputData.ts(end+1,1) = art + crossingIdx(cross)/30000;
                outputData.waveforms(end+1,1:opts.PRE_OFFSET+opts.POST_OFFSET+1) = filterArtifact(art,crossingIdx(cross)-opts.PRE_OFFSET:crossingIdx(cross)+opts.POST_OFFSET);
                outputData.elec(end+1,1) = ch;
            end
            
            
        end
    end

end


function [opts] = configureOptions(optsInput)

    opts.THRESHOLD_MULT = 3.5;
    opts.THRESHOLD_WINDOW_LENGTH = 50;
    opts.MAX_CHAIN_LENGTH = 10;
    opts.PRE_OFFSET = 22;
    opts.POST_OFFSET = 25;
    opts.MAX_AMPLITUDE = 1000;
    opts.NUM_PAD = 400;
    
    [opts.B,opts.A] = butter(6,500/(30000/2),'high');
    
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

