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
    
    spikes_post_stim = cell(size(array_data.numStims));
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
    std_baseline_count = std(all_bin_counts);
    threshold = mean_baseline_count+2*std_baseline_count;
    
    for chan = 1:size(array_data.numStims,1)
        for wave = 1:size(array_data.numStims,2)
            %% collect all unit data into two arrays (one for time, one for stim num)
            spike_times = [];
            stim_nums = [];
            
            for unit = 1:numel(array_data)
                spike_times = [spike_times, array_data.spikeTrialTimes{chan,wave}];
                stim_nums = [stim_nums, array_data.stimData{chan,wave}];
            end
            
            
            %% remove spikes that are not in the window
            post_mask = spike_times > opts.POST_WINDOW(1) & spike_times < opts.POST_WINDOW(2);
            pre_mask = spike_times > opts.PRE_WINDOW(1) & spike_times < opts.PRE_WINDOW(2);
            
            %% count how many spikes for each stimulation
            spikes_post_stim{chan,wave} = zeros(array_data.numStims(chan,wave),1);
            spike_times_post_stim{chan,wave} = array_data.spikeTrialTimes{chan,wave}(post_mask);
            
            spikes_pre_stim{chan,wave} = zeros(array_data.numStims(chan,wave),1);
            
            for stim = 1:array_data.numStims(chan,wave)
                spikes_post_stim{chan,wave}(stim) = sum(stim_nums(post_mask) == stim);
                spikes_pre_stim{chan,wave}(stim) = sum(stim_nums(pre_mask) == stim)/diff(opts.PRE_WINDOW)*diff(opts.POST_WINDOW);
            end
            
            %% compute excitatory statistics based on on (kraskov 2011) or a wilcoxon rank sum test
            
            bins_above_threshold = array_data.binCounts{chan,wave} > threshold & array_data.binEdges{chan,wave}(2:end) > 0 &  array_data.binEdges{chan,wave}(2:end) < 20;
%             is_excitatory(chan,wave) = any(movmean(bins_above_threshold,3) >0.99);
            is_excitatory_p(chan,wave) = signrank(spikes_post_stim{chan,wave},spikes_pre_stim{chan,wave},'tail','right');
        end
    end
        
    % get number of spikes per stim
    num_spikes_post = cellfun(@sum,spikes_post_stim);
    expected_spikes_post = cellfun(@sum,spikes_pre_stim)/diff(opts.PRE_WINDOW)*diff(opts.POST_WINDOW);
    num_stims = cellfun(@numel,spikes_post_stim);
    prob_spikes = (num_spikes_post-expected_spikes_post)./num_stims;

    % determine which conditions to plot based on PW1, PW2, pol. Store amp
    % for those conditions
    amp = zeros(size(prob_spikes));
    pw = zeros(size(prob_spikes,1),2);
    pol = zeros(size(prob_spikes));
    
    keep_mask = zeros(size(prob_spikes)); % determines which conditions to keep
    for cond = 1:numel(array_data.binCounts)
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
            plot(x_data,y_data,'-','marker','.','markersize',16,'color',getColorFromList(1,0))
        end
    else
        %   plot prob_spikes against amp
        plot(amp(keep_mask==1),prob_spikes(keep_mask==1),'-','marker','.','markersize',12,'color',[0.5,0.5,0.5])
    end
    
   output_data.spikes_post_stim = spikes_post_stim;
   output_data.spike_times_post_stim = spike_times_post_stim;
   
   output_data.prob_spike = prob_spikes;
   output_data.keep_mask = keep_mask;
   output_data.amp = amp;
   output_data.pw = pw;
   output_data.pol = pol;
   output_data.is_excitatory_p = is_excitatory_p;
   output_data.mean_baseline_count = mean_baseline_count;
   output_data.std_baseline_count = std_baseline_count;

end