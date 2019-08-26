function [spikes_post_stim,figure_handles] = plotSpikesPerAmplitude(array_data,opts)

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
            spikes_pre_stim{chan,wave} = zeros(array_data.numStims(chan,wave),1);
            
            for stim = 1:array_data.numStims(chan,wave)
                spikes_post_stim{chan,wave}(stim) = sum(stim_nums(post_mask) == stim);
                spikes_pre_stim{chan,wave}(stim) = sum(stim_nums(pre_mask) == stim);
            end
            
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
    keep_mask = zeros(size(prob_spikes)); % determines which conditions to keep
    for cond = 1:numel(array_data.binCounts)
        amp(cond) = array_data.STIM_PARAMETERS(cond).amp1;

        if(array_data.STIM_PARAMETERS(cond).polarity == opts.POL && ...
                array_data.STIM_PARAMETERS(cond).pWidth1 == opts.PW1 && ...
                array_data.STIM_PARAMETERS(cond).pWidth2 == opts.PW2)
            % then keep this condition
            keep_mask(cond) = 1;
        end

    end

    if(opts.MAKE_FIGURE)
        figure_handles = figure();
    end
    
    
    % plot prob_spikes against amp
    plot(amp(keep_mask==1),prob_spikes(keep_mask==1),'-','marker','.','markersize',16)
    

end