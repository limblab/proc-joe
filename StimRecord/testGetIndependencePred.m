%% tests getIndependencePrediction 
% two things
p_list = rand(1000,2);
independence_pred = [];
independence_truth = p_list(:,1)+p_list(:,2) - p_list(:,1).*p_list(:,2);

for i = 1:size(p_list,1)
    independence_pred(i,1) = getIndependencePrediction(p_list(i,:));
end


%% three things
p_list = rand(1000,3);
independence_pred = [];
independence_truth = p_list(:,1)+p_list(:,2)+p_list(:,3) - p_list(:,1).*p_list(:,2)-p_list(:,1).*p_list(:,3)-p_list(:,2).*p_list(:,3) + p_list(:,1).*p_list(:,2).*p_list(:,3);

for i = 1:size(p_list,1)
    independence_pred(i,1) = getIndependencePrediction(p_list(i,:));
end


%% five things
p_list = 0.2*rand(1,5);
num_trials = 100000;
res = zeros(num_trials,5);

for i = 1:numel(p_list)
    res(:,i) = rand(num_trials,1) < p_list(i);
end

sum(sum(res,2)> 0)/num_trials
getIndependencePrediction(p_list)