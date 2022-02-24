function [ mutPop ] = mutatePopulation( population )

global popSize N numPoints numElites numRandom
mutRate = 0.2;

mutPop = population;
r = rand(popSize,1);
numMut = sum(r<mutRate);

for i = 1:numMut
    idxPop = ceil(rand()*(popSize-numElites))+numElites;
    idx = ceil(rand()*N*numPoints);
    str = mutPop(idxPop,:);
    if(str(idx) == '1')
        str(idx) = '0';
    else
        str(idx) = '1';
    end
    mutPop(idxPop,:) = str;
end
% for i = 1:popSize
%     r = rand();
%     if(r < mutRate)
%         idx = ceil(rand()*N*numPoints);
%         str = population{i};
%         if(str(idx) == '1')
%             str(idx) = '0';
%         else
%             str(idx) = '1';
%         end
%         mutPop{i} = str;
%     end
% end
end

