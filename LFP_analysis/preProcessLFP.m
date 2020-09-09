function [ trial_data, lfp_data ] = preProcessLFP(cds, trial_data, input_data )
    % subtract regressed common average from each electrode
    % FFT in 256ms (?) windows
    % calculate power in frequency bands (8-19, 20-69, 70-129, 130-199,
    % 200-300) at each spike bin.... 
    
    t_data = cds.lfp.t;
    lfp_data_pre = cds.lfp{:,2:end};
    
    if(input_data.subtract_common_average)
        common_average = mean(lfp_data_pre,2); % across electrodes
        reg_coeff = common_average\lfp_data_pre;
        
        lfp_data = lfp_data_pre - common_average*reg_coeff;
        
    else
        lfp_data = lfp_data_pre;
    end
    lfp_t = cds.lfp.t;
    lfp_is_baseline = zeros(size(lfp_data,1),1);
    
    % FFT with window length provided, compute power in bands for each
    % channel, store power. Do this for each trial in trial data
    dt = mode(diff(t_data));
    n_points = ceil(input_data.FFT_length/dt/1000);

    for i_trial = 1:numel(trial_data)
        power_data = zeros(size(trial_data(i_trial).pos,1),size(lfp_data,2),size(input_data.power_bands,1));
        cds_trial_idx = find(cds.trials.number == trial_data(i_trial).trial_id);
        
        % find first lfp_idx that corresponds to the current trial
        lfp_idx = find(t_data > cds.trials.startTime(cds_trial_idx)-trial_data(i_trial).bin_size*trial_data(i_trial).idx_startTime,1,'first');
        freqs = 1/dt*(0:(n_points/2))/n_points;
        % loop through each index in the trial, compute power in the
        % corresponding bands and store
        for i_power = 1:size(power_data,1)
            curr_fft = fft(lfp_data(lfp_idx+i_power-n_points:lfp_idx+i_power-1,:));
            curr_power_2_sides = dt/n_points*abs(curr_fft/n_points).^2;
            curr_power_1_side = curr_power_2_sides(1:floor(n_points/2)+1,:);
            curr_power_1_side(2:end-1,:) = curr_power_1_side(2:end-1,:)*2; % double all values except 0 and the last point
            
            for i_band = 1:size(input_data.power_bands,1)
                freq_idx = [find(freqs > input_data.power_bands(i_band,1),1,'first'), ...
                    find(freqs > input_data.power_bands(i_band,2),1,'first')];

                power_data(i_power,:,i_band) = mean(curr_power_1_side(freq_idx(1):freq_idx(2),:));
            end
        end
        
        % store power_data into trial_data
        trial_data(i_trial).lfp_data = power_data;
        
        % get baseline flag based on trial data
        
        
        
    end
    
    

end

