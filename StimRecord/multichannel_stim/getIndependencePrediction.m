function [ independent_p ] = getIndependencePrediction( p )
% array of probabilities for individual events. This function outputs the
% probability of any event occurring if they are all independen
    
    p(p < 0) = 0;
    independent_p = 1-prod((1-p)); % add all together
    
end

