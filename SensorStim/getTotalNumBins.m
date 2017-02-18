function [ totalNumBins ] = getTotalNumBins( common )
% gets the total number of bins in common :D
totalNumBins = 0;
for i = 1:length(common.muscleNames)
    muscleName = common.muscleNames{i};
    totalNumBins = totalNumBins + length(common.(genvarname(muscleName)).counts{1});
end


end

