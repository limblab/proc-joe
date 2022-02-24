function [output_data,figure_handles] = plotSpikesPerCondition(array_data,opts)

    % plots spikes per stimulation for each amplitude condition after subtracting baseline.
    % opts
    % contains PW1, PW2 (length of two phases in us) and POL (polarity of
    % pulse (0 = cathodic-first, 1 = anodic-first). opts should also
    % contain PRE_WINDOW and POST_WINDOW, which are the windows before and
    % after stimulation to look at (in s)
    
    % This function handles a single neuron at a time in case the
    % amplitudes used are different across experiments.
    
    figure_handles = [];
    %% get spikes per stim
    
    num_spikes_post_stim = cell(size(array_data.numStims));
    spikes_pre_stim = cell(size(array_data.numStims));
    spike_times_post_stim = cell(size(array_data.numStims));
    
    is_excitatory_p = zeros(size(array_data.numStims));
    
    % compute mean and std of baseline bin counts
    all_bin_counts = [];
    for chan = 1:size(array_data.numStims,1)
        for wave = 1:size(array_data.numStims,2)
            baseline_idx = [find(array_data.binEdges{chan,wave} > opts.PRE_WINDOW(1)*1000,1,'first'),...
                find(array_data.binEdges{chan,wave} > opts.PRE_WINDOW(2)*1000,1,'first')];
            all_bin_counts = [all_bin_counts, array_data.binCounts{chan,wave}(baseline_idx(1):baseline_idx(2))];
        end
    end
                
    mean_baseline_count = mean(all_bin_counts);
    mean_baseline_fr = mean(all_bin_counts)/(mode(diff(array_data.binEdges{1,1}))/1000);
    std_baseline_count = std(all_bin_counts);
    threshold = mean_baseline_count+2*std_baseline_count;
    post_window_list = [];
    best_window_list = [];
    for chan = 1:size(array_data.numStims,1)
        for wave = 1:size(array_data.numStims,2)
            %% collect all unit data into two arrays (one for time, one for stim num)
            spike_times = [];
            stim_nums = [];
            
            for unit = 1:numel(array_data)
                spike_times = [spike_times, array_data.spikeTrialTimes{chan,wave}];
                stim_nums = [stim_nums, array_data.stimData{chan,wave}];
            end
            
            if(opts.ADJUST_SPIKE_TIMES)
                spike_times = spike_times - opts.WAVEFORM_LENGTH;
            end
            %% remove spikes that are not in the window                
            if(size(opts.POST_WINDOW,1) > 1)
                % find corresponding amp idx and use that window
                post_window_idx = find(array_data.STIM_PARAMETERS(wave).amp2 == opts.AMP_LIST);
                post_mask = spike_times > opts.POST_WINDOW(post_window_idx,1) & spike_times < opts.POST_WINDOW(post_window_idx,2);
                post_window_list(wave,:) = opts.POST_WINDOW(post_window_idx,:);
            else
                post_mask = spike_times > opts.POST_WINDOW(1) & spike_times < opts.POST_WINDOW(2);
                post_window_list = opts.POST_WINDOW(1,:);
            end
            pre_mask = spike_times > opts.PRE_WINDOW(1) & spike_times < opts.PRE_WINDOW(2);
            
            %% count how many spikes for each stimulation
            num_spikes_post_stim{chan,wave} = zeros(array_data.numStims(chan,wave),1);
            spike_times_post_stim{chan,wave} = spike_times(post_mask);
            
            spikes_pre_stim{chan,wave} = zeros(array_data.numStims(chan,wave),1);
            
            for stim = 1:array_data.numStims(chan,wave)
                num_spikes_post_stim{chan,wave}(stim) = sum(stim_nums(post_mask) == stim);
                
                if(size(opts.POST_WINDOW,1) > 1)
                    post_window_diff = diff(opts.POST_WINDOW(post_window_idx,:)); % found in first if statement (remove spikes block)
                else
                    post_window_diff = diff(opts.POST_WINDOW(1,:));
                end
                spikes_pre_stim{chan,wave}(stim) = sum(stim_nums(pre_mask) == stim);
            end
            
            if(opts.MAXIMIZE_RESPONSE_PROB) % find window with most spikes instead of using spike_times(post_mask)
                best_response_prob = -1000;
                best_window = [-1,-1];
                
                for window_size = 0.5/1000:0.5/1000:diff(opts.POST_WINDOW)
                    for start_time = opts.POST_WINDOW(1):window_size:opts.POST_WINDOW(2)
                        window_test = [start_time,start_time+window_size];
                        num_spikes = numel(spike_times(spike_times > window_test(1) & spike_times < window_test(2)));
                        response_prob = (num_spikes-sum(spikes_pre_stim{chan,wave})/diff(opts.PRE_WINDOW)*diff(window_test))/array_data.numStims(chan,wave);
                        
                        if(response_prob > best_response_prob)
                            best_window = window_test;
                            best_response_prob = response_prob;
                        end
                    end
                end
                
                % override parameters
                post_mask = spike_times > best_window(1) & spike_times < best_window(2);
                spike_times_post_stim{chan,wave} = spike_times(post_mask);
                num_spikes_post_stim{chan,wave} = zeros(array_data.numStims(chan,wave),1);
                for stim = 1:array_data.numStims(chan,wave)
                    num_spikes_post_stim{chan,wave}(stim) = sum(stim_nums(post_mask) == stim);
                end
                
                best_window_list(wave,:) = best_window;
            end
            %% compute excitatory statistics based on on (kraskov 2011) or a wilcoxon rank sum test
            
            bins_above_threshold = array_data.binCounts{chan,wave} > threshold & array_data.binEdges{chan,wave}(2:end) > 0 &  array_data.binEdges{chan,wave}(2:end) < 20;
%             is_excitatory(chan,wave) = any(movmean(bins_above_threshold,3) >0.99);
            is_excitatory_p(chan,wave) = signrank(num_spikes_post_stim{chan,wave},spikes_pre_stim{chan,wave}/diff(opts.PRE_WINDOW).*diff(post_window_list'),'tail','right');
        end
    end
        
    if(opts.MAXIMIZE_RESPONSE_PROB)
        post_window_list = best_window_list;
    end
    % get number of spikes per stim
    num_spikes_post = cellfun(@sum,num_spikes_post_stim);
    expected_spikes_post = cellfun(@sum,spikes_pre_stim)/diff(opts.PRE_WINDOW).*diff(post_window_list');
    num_stims = cellfun(@numel,num_spikes_post_stim);
    prob_spikes = (num_spikes_post-expected_spikes_post)./num_stims;
%     prob_spikes = num_spikes_post./num_stims;
    
    % determine which conditions to plot based on PW1, PW2, pol. Store amp
    % for those conditions
    amp = zeros(size(prob_spikes));
    pw = zeros(size(prob_spikes,1),2);
    pol = zeros(size(prob_spikes));
    
    keep_mask = zeros(size(prob_spikes)); % determines which conditions to keep
    for cond = 1:size(array_data.binCounts,2)
        amp(cond) = array_data.STIM_PARAMETERS(cond).amp2;
        pw(cond,:) = [array_data.STIM_PARAMETERS(cond).pWidth1,array_data.STIM_PARAMETERS(cond).pWidth2];
        pol(cond) = array_data.STIM_PARAMETERS(cond).polarity;
        
        if((ismember(amp(cond),opts.AMP) || isempty(opts.AMP)) && ...
            (ismember(pol(cond),opts.POL) || isempty(opts.POL)) && ...
            (ismember(pw(cond,1),opts.PW1) || isempty(opts.PW1)) && ...
            (ismember(pw(cond,2),opts.PW2) || isempty(opts.PW2)))
                
            % then keep this condition
            keep_mask(cond) = 1;
        end

    end

    if(opts.MAKE_FIGURE)
        figure_handles = figure();
    end
    
    if(opts.PLOT_PULSE_WIDTH)
        x_data = [1:1:sum(keep_mask)];
        y_data = prob_spikes(keep_mask == 1);
        if(numel(y_data) == numel(x_data)) % some data in Han does not have the pulse width conditions
            plot(x_data,y_data,'-','marker','none','markersize',16,'color',getColorFromList(1,0),'linewidth',1)
        end
    else
        %   plot prob_spikes against amp
        
        plot(amp(keep_mask==1),prob_spikes(keep_mask==1),'-','marker','.','markersize',12,'color',opts.COLOR,'linewidth',1)
    end
    
   output_data.num_spikes_post_stim = num_spikes_post_stim;
   output_data.spike_times_post_stim = spike_times_post_stim;
   
   output_data.prob_spike = prob_spikes;
   output_data.keep_mask = keep_mask;
   output_data.amp = amp;
   output_data.pw = pw;
   output_data.pol = pol;
   output_data.is_excitatory_p = is_excitatory_p;
   output_data.mean_baseline_count = mean_baseline_count;
   output_data.mean_baseline_fr = mean_baseline_fr;
   output_data.std_baseline_count = std_baseline_count;
   
   output_data.post_window_list = post_window_list;

end