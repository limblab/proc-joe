function [ common ] = getAllNeuronCounts( common )
% populates the neuronCounts matrix in common by counting all neurons in
% common.<muscle name>.<counts>. Also creates and populates
% common.muscleIdx which contains the starting, ending, and zero index for each
% muscle in common.neuronCounts ordered by common.muscleNames.

common.muscleIdx = cell(length(common.muscleNames),4);
for i = 1:length(common.muscleNames) % for each muscle
    
   % get start idx, end idx, and zero idx
   if(i == 1)
       common.muscleIdx{i,1} = common.muscleNames{i};
       common.muscleIdx{i,2} = 1; % start idx
       common.muscleIdx{i,4} = common.muscleIdx{i,2} + ...
           length(common.(genvarname(common.muscleNames{i})).counts{1})-1; % end idx
       zeroEdge = find(common.(genvarname(common.muscleNames{i})).edges{i}==0);
       common.muscleIdx{i,3} = common.muscleIdx{i,2} + zeroEdge - 0.5 - 1; % zero idx
   else
       common.muscleIdx{i,1} = common.muscleNames{i};
       common.muscleIdx{i,2} = common.muscleIdx{i-1,4}+1; % start idx
       common.muscleIdx{i,4} = common.muscleIdx{i,2} + ...
           length(common.(genvarname(common.muscleNames{i})).counts{1})-1; % end idx
       zeroEdge = find(common.(genvarname(common.muscleNames{i})).edges{i}==0);
       common.muscleIdx{i,3} = common.muscleIdx{i,2} + zeroEdge - 0.5 - 1; % zero idx
   end
   
   % populate common.neuronCounts for each neuron
   startIdx = common.muscleIdx{i,2};
   endIdx = common.muscleIdx{i,4};
   
   for j = 1:length(common.muscles)
        common.neuronCounts(j,startIdx:endIdx) = common.(genvarname(common.muscleNames{i})).counts{j};
   end
    
end


end

