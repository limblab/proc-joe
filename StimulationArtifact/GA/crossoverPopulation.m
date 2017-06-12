function [ crossPop ] = crossoverPopulation( population )

global popSize N numPoints numElites 

pc_rate = 0.8;
crossPop = population;
count = 1;
r = rand(popSize-numElites,1);
numCross = sum(r<pc_rate);
for i = 1:numCross
    idx = ceil(rand()*(popSize-numElites)) + numElites;
    idxOther = idx;
    while(idxOther==idx)
        idxOther = ceil(rand()*(popSize-numElites))+numElites;
    end
    splitSpot = ceil(rand()*(N*numPoints-4))+2; 
    str1 = population(idx,:);
    str2 = population(idxOther,:);
    crossPop(idxOther,:) = strcat(str1(1:splitSpot),str2(splitSpot+1:end));
    crossPop(idx,:) = strcat(str2(1:splitSpot),str1(splitSpot+1:end));
end

% for i = 1:popSize
%     r = rand();
%     if(r < pc_rate)
%         %% cross this and a random (not this one) chromosome
%         % replace both
%         idxOther = i;
%         while(idxOther==i)
%             idxOther = ceil(rand()*popSize);
%         end
%         splitSpot = ceil(rand()*(N*numPoints-2))+1;
%         str1 = population{i};
%         str2 = population{idxOther};
%         crossPop{idxOther} = strcat(str1(1:splitSpot),str2(splitSpot+1:end));
%         crossPop{i} = strcat(str2(1:splitSpot),str1(splitSpot+1:end));
%     end
% end

end

