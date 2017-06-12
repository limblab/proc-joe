function [ts, waves, chan, numCrossings, shouldFind] = getThresholdCrossingData(data, thresholdMult, preOffset,postOffset,idxPlaced,art)
lengthWave = postOffset + preOffset + 1; % plus one for the 0 idx
ts = zeros(1000,1);
waves = zeros(1000,lengthWave);
chan = zeros(1000,1);
numCrossings = 1;
shouldFind = zeros(96,1);

for stim = 1:size(data.artifact,2)
    for elec = 1:size(data.artifact,1)
%         elec = 1;
%         stim=stim+1;
        stimData = squeeze(data.artifact(elec,stim,1:end));
        threshold = abs(rms(stimData(end-10:end))*thresholdMult);
        thresholdCrossings = find(stimData(:)>threshold);
%         figure;
%         plot(stimData)
%         hold on
%         plot([0 280],[threshold threshold])
%         plot([0 280],[-threshold -threshold])

        % check for correct index 
%         flag = 0;
%         for cross = 1:numel(thresholdCrossings)
%             if(thresholdCrossings(cross) > idxPlaced-2 && thresholdCrossings(cross) < idxPlaced+2)
%                 flag = 1;
%             end
%         end
%         shouldFind(elec) = shouldFind(elec) + flag;

        
        % remove chains -- find min in chain and keep that        
        idx = 2;
        chain = [1];
        crossingsKeep = [];
        while idx <= numel(thresholdCrossings)
            if(thresholdCrossings(idx) == thresholdCrossings(idx-1)+1) % store in chain
                chain = [chain;idx];
            elseif(~isempty(chain)) % broke a chain, store minidx, update idx, empty chain
                [~,maxIdx] = max(stimData(thresholdCrossings(chain)));
                if(isempty(crossingsKeep))
                    crossingsKeep = [thresholdCrossings(maxIdx+chain(1)-1)];
                else
                    crossingsKeep = [crossingsKeep;thresholdCrossings(maxIdx+chain(1)-1)];
                end
                chain = [idx];
            end
            idx = idx+1;
        end
        if(numel(thresholdCrossings) > 0)
            thresholdCrossings = [crossingsKeep;thresholdCrossings(end)];
        end
        % remove potential artifacts and too close to beginning
        crossingsMask = ones(numel(thresholdCrossings),1);
        for cross = 1:numel(thresholdCrossings)
            if(stimData(thresholdCrossings(cross)) > 1000 || ...
                    ~(thresholdCrossings(cross)+postOffset <= numel(stimData) && thresholdCrossings(cross)-preOffset > 0))
                crossingsMask(cross) = 0;
            end
        end
        thresholdCrossings = thresholdCrossings(crossingsMask(:) == 1);

        % go through and weed out ones that are too close - prioritize
        % backwards in time
        crossingsMask = ones(numel(thresholdCrossings),1);
        for cross = numel(thresholdCrossings):-1:2
            if(crossingsMask(cross) == 1) % check time beforehand to see if one is too close
                crossCheck = cross-1;
                while crossCheck >= 1 && thresholdCrossings(crossCheck) >= thresholdCrossings(cross) - max(preOffset,postOffset)
                    crossingsMask(crossCheck) = 0;
                    crossCheck = crossCheck-1;
                end
            end
        end  
        thresholdCrossings = thresholdCrossings(crossingsMask(:) == 1);
        
        % store thresholdCrossing data
        for cross = 1:numel(thresholdCrossings)
            ts(numCrossings) = 1 + (art-1)*10 + ((stim-1)*0.01) + thresholdCrossings(cross)/30000; % this is in seconds
            chan(numCrossings) = elec;
            waves(numCrossings,:) = stimData(thresholdCrossings(cross)-preOffset:thresholdCrossings(cross)+postOffset);
            numCrossings = numCrossings + 1;

            if(numCrossings > 0.67*numel(ts))
                ts = [ts;zeros(1000,1)];
                waves = [waves; zeros(1000,lengthWave)];
                chan = [chan; zeros(1000,1)];
            end
        end
    end 
end

end


% old thresholding method

% % move backwards in time and check for artifacts and chains
%         cross = numel(thresholdCrossings);
%         while cross >= 1    
%             % determine if in chain, if so move cross to min in chain,
%             % remove all others
%             chain = [cross];
%             while chain(1) >= 2 && thresholdCrossings(chain(1)-1) + 1 == thresholdCrossings(chain(1))
%                 chain = [chain(1)-1;chain(:)];
%             end
%             [~,crossUse] = min(stimData(thresholdCrossings(chain)));
%             crossUse = crossUse+chain(1)-1;
%             if(thresholdCrossings(crossUse)+postOffset <= numel(stimData) && thresholdCrossings(crossUse)-preOffset > 0)
%                 ts(numCrossings) = ((stim-1)*0.01) + thresholdCrossings(crossUse)/30000; % this is in seconds
%                 chan(numCrossings) = elec;
%                 waves(numCrossings,:) = stimData(thresholdCrossings(crossUse)-preOffset:thresholdCrossings(crossUse)+postOffset);
%                 numCrossings = numCrossings + 1;
% 
%                 if(numCrossings > 0.67*numel(ts))
%                     ts = [ts;zeros(1000,1)];
%                     waves = [waves; zeros(1000,lengthWave)];
%                     chan = [chan; zeros(1000,1)];
%                 end
%                 % update cross -- move backwards enough in time
%                 newCross = cross;
%                 while newCross >= 1 && thresholdCrossings(crossUse)-thresholdCrossings(newCross)<max(preOffset,postOffset)
%                     newCross = newCross - 1;
%                 end
%                 cross = newCross;
%             else
%                 cross = cross - 1;
%             end
%         end