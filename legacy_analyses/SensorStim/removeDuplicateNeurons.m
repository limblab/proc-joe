function [ out_common ] = removeDuplicateNeurons( common )
% removes duplicate neurons from common by not porting them over to
% out_common. Checks if channel, id, electrode match for a neuron in common
% to any neuron in out_common. If match, no move over, else move over. 

out_common.muscles = {};
out_common.channels = [];
out_common.IDs = [];
out_common.electrodes = {};
out_common.neuronNumbers = [];

for i = 1:length(common.muscles) % for each neuron in common
    moveOver = 1;
    if(length(out_common.muscles) > 0)
        for j = 1:length(out_common.muscles) % for each neuron in out_common
            if(common.channels(i,1) == out_common.channels(j,1) && common.IDs(i,1) == out_common.IDs(j,1) && ... 
                    strcmp(common.electrodes{i,1},out_common.electrodes{j,1}))
                moveOver = 0;
            end
        end
    end
    if(moveOver)
        out_common.muscles{end+1,1} = common.muscles{i};
        out_common.electrodes{end+1,1} = common.electrodes{i};
        out_common.channels(end+1,1) = common.channels(i);
        out_common.IDs(end+1,1) = common.IDs(i);
        out_common.neuronNumbers(end+1,1) = common.neuronNumbers(i);
    end
end

end

