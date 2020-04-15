function [output_data] = bootstrapPsychData(psych_data,input_data)

    % bootstrapping on psych data
    % get fit for bump data and all stim data again
    output_data = [];
    
    % determine how many trials total, then randomize idx's and which entry
    % they refer to
    total_trials = 0;
    psych_data_idx = []; % stores condition in psych_data -- psych_data{i}(j) is converted to 
    counter = 1;
    for i = 1:numel(psych_data)
        for j = 1:numel(psych_data{i})
            total_trials = total_trials + numel(psych_data{i}(j).trial_ids);
            psych_data_idx = [psych_data_idx, counter + zeros(1,numel(psych_data{i}(j).trial_ids))];
            counter = counter+1;
        end
    end
    
    trial_idx = repmat(1:total_trials,input_data.num_bootstrap,1);
    psych_data_idx = repmat(psych_data_idx,input_data.num_bootstrap,1);
    % sample trial_idx with replacement num_bootstrap times
    for boot = 1:size(trial_idx,1)
        [trial_idx(boot,:),sample_idx] = datasample(trial_idx(boot,:),size(trial_idx,2),'Replace',true);
        psych_data_idx(boot,:) = psych_data_idx(boot,sample_idx);
    end
    
    % apply sampling to psych_data_idx
    
    % for each condition, compute bump_correct, bump_total, and fit psych
    % curve as well as get point of subjective equality
    
    
    for boot = 1:input_data.num_bootstrap
        % sample trial_idx with replacement num_bootstrap times
        
        
        
    end

end