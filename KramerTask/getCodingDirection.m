function [output_data] = getCodingDirection(data,input_data)

    % function returns a num neurons x time (however many bins are given)
    % vector of the coding direction -- defined as r_right - r_left for
    % each neuron
    output_data = {};
    
    training_trial_idx = datasample(1:1:size(data.spike_data,1),floor(input_data.sample_rate*size(data.spike_data,1)),'Replace',false);
    training_trial_mask = zeros(size(data.spike_data,1),1);
    training_trial_mask(training_trial_idx,1) = 1;
    
    cd = squeeze(mean(data.spike_data(data.reached_0_deg == 1 & training_trial_mask == 1,:,:),1) - ...
        mean(data.spike_data(data.reached_0_deg == 0 & training_trial_mask == 1,:,:),1));

    output_data.cd = cd;
    output_data.mean_cd = mean(cd,2);
    output_data.training_trial_mask = training_trial_mask;
end