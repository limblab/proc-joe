function [output_data] = getCodingDirection(data,input_data)

    % function returns a num neurons x time (however many bins are given)
    % vector of the coding direction -- defined as r_right - r_left for
    % each neuron
    output_data = {};
    
    % only use bump directions within a specified range
    if(isfield(input_data,'bump_range'))
        bump_mask = data.bump_dirs > input_data.bump_range(1) & data.bump_dirs < input_data.bump_range(2);
    else
        bump_mask = ones(size(data.bump_dirs));
    end
    
    valid_trials = find(bump_mask);
    
    % split data up based on training trials
    training_trial_idx = datasample(1:1:size(valid_trials,1),floor(input_data.sample_rate*size(valid_trials,1)),'Replace',false);
    training_trial_mask = zeros(size(data.spike_data,1),1);
    training_trial_mask(valid_trials(training_trial_idx),1) = 1;
    
    
    cd = squeeze(mean(data.spike_data(data.reached_0_deg == 1 & training_trial_mask == 1 & bump_mask == 1,:,:),1) - ...
        mean(data.spike_data(data.reached_0_deg == 0 & training_trial_mask == 1 & bump_mask == 1,:,:),1));

    output_data.cd = cd;
    output_data.mean_cd = mean(cd,2);
    output_data.training_trial_mask = training_trial_mask;
    output_data.bump_mask = bump_mask;
    output_data.bump_range = input_data.bump_range;
end