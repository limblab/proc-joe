function [output_data] = makeStimArrayWrapper(pattern_data,stim_params)
    % given an IPI and a set of desired frequencies (pattern data),
    % generate a stim array for each pattern
    % stim_params contains IPI, max_freq, train_length

    for i_pattern = 1:numel(pattern_data.pattern)
        desired_freqs = pattern_data.pattern{i_pattern}.stim_norm * stim_params.max_freq;
        
        [stim_array, metrics] = makeStimArray(desired_freqs,stim_params);
    end

    output_data = [];


end




function [stim_array,metrics] = makeStimArray(freqs,stim_params)
    % freqs in Hz
    % stim_params contains IPI, max_freq, train_length
    stim_array = []; metrics =[];
    % define useful variables
    best_score = -1;
    best_stim_array = 0;
    
    % initialize stim array that satisfy freqs
    stim_array = zeros(numel(freqs), ceil(stim_params.train_length/(stim_params.IPI)));
    
    % for a bunch of generations
        
        % score each set
    
        % if one set satisfies the requirements, end and output that set
        
        % else, randomly select best stim_arrays and perturb
        
        % repeat for next generation


end


function [score] = getStimArrayScore(desired_freqs, stim_array)

    % score is defined as
    % sum of inter-pulse interval minus desired inter-pulse interval
    % 1000000 if more than 16 electrodes are stimmed in a single time step
    score = 0;
    for i_elec = 1:size(stim_array,1)
        IPI_list = diff(find(stim_array(i_elec,:)));
        
        desired_IPI = 1000/desired_freqs(i_elec);
        score = score + sum(abs(IPI_list-desired_IPI));
    end
    % add if more than 16 elecs are stimmed in a single time step
    score = score + 1000000*any(sum(stim_array,1) > 16);

end