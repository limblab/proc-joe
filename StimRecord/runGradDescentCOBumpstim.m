function [ guess,count ] = runGradDescent( neuronMeanWaveFilt,b,a,varargin )
%% line search currently does not work
lineSearch = 0;
beta = 0;
alpha = 1;
for i = 1:2:numel(varargin)
    switch varargin{i}
        case 'lineSearch'
            lineSearch = 1;
            temp = varargin{i+1};
            alpha = temp(1);
            beta = temp(2);
    end
end

count = 1;
maxCount = 5000;
threshold = 0.0000001;
alpha = 100;
step = 1;

guess = zeros(maxCount,numel(neuronMeanWaveFilt));
guess(count,:) = rand(1,numel(neuronMeanWaveFilt))*800-400;

gradJ = threshold + 1;
while(count < maxCount && norm(gradJ,2) > threshold)
    % gradient
    gradJ = (guess(count,:) + alpha*(fliplr(filter(b,a,fliplr(guess(count,:))))-neuronMeanWaveFilt)*b(1)/a(1))/(alpha+1);
    if(lineSearch)
        n = 0;
        gamma = beta^n;
        while(norm((fliplr(filter(b,a,fliplr(guess(count,:)-gamma*gradJ)))-neuronMeanWaveFilt)) > ... 
                norm((fliplr(filter(b,a,fliplr(guess(count,:))))-neuronMeanWaveFilt) - alpha*gamma*gradJ))
            n = n+1;
            gamma = beta^n;
        end
        guess(count+1,:) = guess(count,:)- gamma*gradJ;
    else
        guess(count+1,:) = guess(count,:) - gradJ*step;
    end
    
    count = count+1;
end

end

