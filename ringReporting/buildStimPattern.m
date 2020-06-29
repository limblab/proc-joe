function [ pattern_data ] = buildStimPattern( neural_data,chan_list, input_data )
    % this function generates stimulation patterns based on input_data
    % input_data contains
    %   constant_scramble;
    %   num_scrambles;
    %   min_msd_percentile
    %   use_pd = will use cosine instead of treat neural data as a z_score
    %   neural_data_idx_list = list of indices in neural_data. Will make
    %       twice this number of patterns (1 bio, 1 nonbio)
    %   directions = cell array of list of directions. first entry
    %       corresponds to first entry in neural_data_idx_list, and so on
    %   num_patterns = size(neural_data_idx_list,1)
    %   num_chans , if empty or <=0, will use all rows in neural_data
    %   max_data
    %   min_data , maps data between these 0 and 1 using these two
    
    % set up pattern data struct
    pattern_data = [];
    pattern_data.num_patterns = numel(tgt_dirs)*2; % must be even (pairs of bio and nonbio)
        
    pattern_data.pred_dirs = [];
    pattern_data.is_biomimetic = [];
    pattern_data.pattern = [];

    % make biomimetic patterns 
    for i_pattern = 1:input_data.num_patterns
        if(~input_data.use_pd) % can do both tgt_dir and bump_dir simultaneously
            % sample channels
            chans_use = chan_list(abs(neural_data(chan_mask,td_idx,bd_idx)) < 1000 & ~isnan(z_score(chan_mask,td_idx,bd_idx)));
            chans = datasample(chans_use,pattern_data.num_chans,'Replace',false); % actual channel numbers
            
            
        else % currently only setup for a single direction
            % sample channels

            % get PD projection for channel list
            PD_proj = cos(input_data.directions{1}(neural_data_idx_list(i_pattern,1) - neural_data);
            PD_proj(PD_proj < 0) = 0;
            
            pattern_data.pred_dirs(end+1) = tgt_dirs(i_tgt);
            bio_pattern = []; 
            bio_pattern.stim_norm = PD_proj; 
            bio_pattern.chans = chan_list(sample_data.groups(:,group_idx),1); % convert idx into channel number
            pattern_data.pattern{end+1} = bio_pattern;
            pattern_data.is_biomimetic(end+1) = 1;
        end
    end
    
end

