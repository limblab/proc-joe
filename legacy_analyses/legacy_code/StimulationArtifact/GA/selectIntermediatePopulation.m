function [ intPop ] = selectIntermediatePopulation( population, allFit )
% return an intermediate population where each member is chosen from the
% previous population with probability proportional to its fitness (in this
% inverse fitness)

global popSize numRandom N numPoints

intPop = population;

for i = 1:popSize
    idx1 = ceil(rand()*popSize);
    idx2 = idx1;
    while(idx2 == idx1)
        idx2 = ceil(rand()*popSize);
    end
    better = idx2;
    worse = idx1;
    if(allFit(idx1,1) < allFit(idx2,1))
        better = idx1;
        worse = idx2;
    end
    r = rand();
    if(r < 0.9)
        intPop(i,:) = population(better,:);
    else
        intPop(i,:) = population(worse,:);
    end
end

[~, bestIdx] = min(allFit);
intPop(1,:) = population(bestIdx,:);
intPop(end-numRandom+1:end,:) = num2str(round(rand(numRandom,N*numPoints)),'%d');

end

