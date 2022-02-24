function [ arrayLabels ] = labelArray( common,mapfile )
% attempts to label the pins on the array based on the neurons in common
% and the given mapfile

arrayMap = loadMapFile(mapfile);
% col - 0 based column from left to right					
% row - 0 based row from bottom to top	

arrayLabels = cell(10);
for i = 1:length(common.neuronLabels)
    for j = 1:size(arrayMap,1)
        if(arrayMap.chan(j) == common.channels(i))
            arrayLabels{arrayMap.row(j),arrayMap.col(j)} = common.neuronLabels{i};
        end
    end
end

end

