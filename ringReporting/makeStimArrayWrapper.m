function [output_data] = makeStimArrayWrapper(pattern_data,stim_params)
    % given an IPI and a set of desired frequencies (pattern data),
    % generate a stim array for each pattern
    % stim_params contains IPI, max_freq, train_length
    output_data = [];
    output_data.stim_array = cell(numel(pattern_data.pattern),1);
    output_data.chans = cell(numel(pattern_data.pattern),1);
    output_data.metrics = cell(numel(pattern_data.pattern),1);
    
    for i_pattern = 1:numel(pattern_data.pattern)
        desired_freqs = pattern_data.pattern{i_pattern}.stim_norm*(stim_params.max_freq-stim_params.min_freq) + stim_params.min_freq;
        chans = pattern_data.pattern{i_pattern}.chans;
        keep_mask = desired_freqs > 0;
        
        desired_freqs = desired_freqs(keep_mask);
        chans = chans(keep_mask);
        [stim_array] = makeStimArray(desired_freqs,stim_params);
        
        output_data.stim_array{i_pattern} = stim_array;
        output_data.chans{i_pattern} = chans;
    end



end




function [best_stim_array] = makeStimArray(freqs,stim_params)
    % freqs in Hz
    % stim_params contains IPI, max_freq, train_length
    % define useful variables
    best_loss = -1;
    best_stim_array = 0;
    num_gens = 2000;
    num_stim_arrays = 100; % make it even please
    
    % initialize stim array that satisfy freqs
    stim_array = zeros(num_stim_arrays, numel(freqs), ceil(stim_params.train_length/(stim_params.IPI)));
    for i_elec = 1:size(stim_array,2)
        start_idx = ceil(i_elec/16);
        stim_steps = round(start_idx:1/freqs(i_elec)/stim_params.IPI:size(stim_array,3));
        stim_array(:,i_elec,stim_steps) = 1;
    end
    
    % for a bunch of generations
    for i_gen = 1:num_gens    
        % score stim_array
        loss_list = getStimArrayScore(freqs, stim_array);
        
        [min_loss,min_loss_idx] = min(loss_list);
        % store score if it's best
        if(min_loss < best_loss || best_loss < 0)
            best_loss = min_loss;
            best_stim_array = squeeze(stim_array(min_loss_idx,:,:));
        end
        
        % randomly perturb stim_array
        stim_array = perturbStimArray(stim_array,loss_list);
        
    end % repeat for next generation
    
    
    % get metrics based on best_stim_array
    
    
end


function [loss_list] = getStimArrayScore(desired_freqs, stim_array)

    % score is defined as
    % sum of inter-pulse interval minus desired inter-pulse interval
    % 1000000 if more than 16 electrodes are stimmed in a single time step
    loss_list = zeros(size(stim_array,1),1);
    
    for i_array = 1:size(stim_array,1)
        for i_elec = 1:size(stim_array,2)
            IPI_list = diff(find(stim_array(i_array,i_elec,:)));
            desired_IPI = 1000/desired_freqs(i_elec);
            loss_list(i_array) = loss_list(i_array) + sum(abs(IPI_list-desired_IPI));
        end
        % add if more than 16 elecs are stimmed in a single time step
        loss_list(i_array) = loss_list(i_array) + 1000000*any(sum(stim_array(i_array,:,:),2) > 16);    
    end
    

end


function [stim_array] = perturbStimArray(stim_array,loss_list)

    % perturb all
    pert_types = rand(size(stim_array,1),1);
    for i_array = 1:size(stim_array,1)
        elec_idx = ceil(rand()*size(stim_array,2));
        if(pert_types(i_array) < 0.33) % swap a '1' and '0'
            one_idx = find(stim_array(i_array,elec_idx,:));
            one_idx = one_idx(ceil(rand()*numel(one_idx)));
            zero_idx = find(stim_array(i_array,elec_idx,:) == 0);
            zero_idx = zero_idx(ceil(rand()*numel(zero_idx)));
            stim_array(i_array,elec_idx,one_idx) = 0;
            stim_array(i_array,elec_idx,zero_idx) = 1;            
        else % shift entire row by 1
            dir = round(rand())*2 - 1;
            stim_array(i_array,elec_idx,:) = circshift(stim_array(i_array,elec_idx,:),dir);
        end
        
    end

end