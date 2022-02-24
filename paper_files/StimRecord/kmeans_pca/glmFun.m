%% glm fun

%% need to make an array of inputs
input_data = [];


%% need to make an array of outputs
output_data = [];
% spike rates, binned per trial
bin_size = 1; % ms
bE = (arrayData{1}.bE{1}(1):bin_size:arrayData{1}.bE{1}(end));
for arr = 1:numel(arrayData)
    for ampIdx = 1:numel(arrayData{arr}.spikeTrialTimes)
        [bC,bE] = histcounts(arrayData{arr}.spikeTrialTimes{ampIdx}*1000,bE);
        output_data = [output_data; bC'];
    end
end
%% then run glmfit

[b,dev,stats] = glmfit(input_data,output_data,'poisson');