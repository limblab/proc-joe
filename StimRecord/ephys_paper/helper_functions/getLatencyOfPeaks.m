function [output_data] = getLatencyOfPeaks(array_data,input_data)


    peak_size = [];
    peak_num = [];
    lat_list = [];
    
    std_list = [];
    width_list = [];
    unit_idx = [];
    amp_list = [];
    
    num_peaks_amp = [];
    num_peaks_amp_assume = [];
    num_peaks_missed = zeros(numel(array_data),1);
    
    delta_lat_list = [];
    delta_amp_list = [];
    delta_unit_idx = [];
    
    % if it's the model, store diam_idx, cell_idx, clone num
    diam_list = [];
    cell_list = [];
    clone_list = [];
    
    % get bin_edges from input_data
    if(input_data.use_gauss_filter == 0)
        bin_edges = input_data.peak_window(1):input_data.bin_size:input_data.peak_window(2);
    else
    
        % generate kernel if using gauss filter later
        % kernel half length is 3·SD out
        kernel_hl = ceil( 3 * input_data.kernel_SD / (input_data.dt) );
        % create the kernel --it will have length 2*kernel_hl+1
        kernel = normpdf( -kernel_hl*(input_data.dt) : ...
            input_data.dt : kernel_hl*(input_data.dt), ...
            0, input_data.kernel_SD );
        % compute normalization factor --this factor depends on the number of taps
        % actually used
        bin_edges = input_data.peak_window(1):input_data.dt:input_data.peak_window(2);
        nm = conv(kernel,ones(1,numel(bin_edges)-1)');
    end
    bin_centers = bin_edges(1:end-1) + mode(diff(bin_edges)/2);
    
    for i_unit = 1:numel(array_data)
        for i_amp = 1:numel(input_data.amp_list)
            amp_idx = find(input_data.amp_list(i_amp) == [array_data{i_unit}.STIM_PARAMETERS.amp1],1,'first');
            
            if(~isempty(amp_idx))
                
%                     spike_mask = array_data{i_unit}.spikeTrialTimes{amp_idx} >= input_data.windows(i_window,1) & ...
%                         array_data{i_unit}.spikeTrialTimes{amp_idx} <= input_data.windows(i_window,2);
                
                spike_mask = array_data{i_unit}.spikeTrialTimes{amp_idx} >= input_data.peak_window(1) & ...
                        array_data{i_unit}.spikeTrialTimes{amp_idx} <= input_data.peak_window(2);  
                    
                spike_times = array_data{i_unit}.spikeTrialTimes{amp_idx}(spike_mask);

                % find peaks
                if(input_data.use_gauss_filter)
                    spike_times_use = conv(kernel,histcounts(spike_times,bin_edges))./nm';
                    spike_times_use = spike_times_use(kernel_hl+1:end-kernel_hl);
                else
                    spike_times_use = histcounts(spike_times,bin_edges);
                end
                
                [pks,locs,widths,proms] = findpeaks(spike_times_use,... % data,
                    'MinPeakDistance',7,'MinPeakProminence',1.0,'NPeaks',4); 
                
                % get time window to count spikes in 
                % do this by drawing a line from peak down based on half
                % width height
                
                stds_temp = nan(1,numel(pks));
                for i_peak = 1:numel(pks)
                    half_prom = proms(i_peak)/2;
                    
                    % find first point below half_prom starting at peak and
                    % moving to the right
                    right_bound = find(spike_times_use(locs(i_peak):end) < half_prom,1,'first') + locs(i_peak) - 1;
                    
                    % find first point below half_prom starting at peak and
                    % moving to the left
                    left_bound = find(spike_times_use(1:locs) < half_prom,1,'last');
                    
                    % extend lines to 0 (line from peak through left or
                    % right bound) to get window to count spikes
                    window = (locs(i_peak) + [-(locs(i_peak)-left_bound)*pks(i_peak)/(pks(i_peak)-spike_times_use(left_bound)),...
                        -(locs(i_peak)-right_bound)*pks(i_peak)/(pks(i_peak)-spike_times_use(right_bound))])*input_data.dt;
                    
                    
                    % get std of spike times in this window
                    stds_temp(i_peak) = std(spike_times(spike_times >= window(1) & spike_times <= window(2)));
                end
                
                
                
                
                % store peak location, width, and size
                peak_size(end+1:end+numel(pks)) = pks;
                peak_num(end+1:end+numel(pks)) = 1:1:numel(pks);
                lat_list(end+1:end+numel(pks)) = bin_centers(locs);
                width_list(end+1:end+numel(pks)) = widths*input_data.dt;
                std_list(end+1:end+numel(pks)) = stds_temp;
                
                unit_idx(end+1:end+numel(pks)) = i_unit;
                amp_list(end+1:end+numel(pks)) = input_data.amp_list(i_amp);
                
                num_peaks_amp(end+1,:) = [numel(pks),input_data.amp_list(i_amp)];
                % if it's the model, store diam_idx, cell_idx, clone num
                if(input_data.is_model)
                    diam_list(end+1:end+numel(pks)) = array_data{i_unit}.diam;
                    cell_list(end+1:end+numel(pks)) = array_data{i_unit}.cell_id;
                    clone_list(end+1:end+numel(pks)) = array_data{i_unit}.clone_num;
                end

                % if experiment, check to see if 100uA condition misses an assumed
                % peak at ~1ms due to the artifact. Add result to an
                % assumed_variable.
                if(~input_data.is_model && input_data.amp_list(i_amp)==50)
                    % store 50uA data
                    first_50_peak = min(bin_centers(locs));
                end
                if(~input_data.is_model && input_data.amp_list(i_amp)==100)                    
                    % if 50uA has a peak and 100uA does not, assume it is
                    % there. 
                    first_100_peak = min(bin_centers(locs));
                    if(isempty(first_50_peak))
                        num_peaks_missed(i_unit) = 0;
                    elseif(isempty(first_100_peak))
                        num_peaks_missed(i_unit) = 1;
                    elseif(first_100_peak > 0.0019 && first_50_peak < 0.0014)
                        num_peaks_missed(i_unit) = 1;
                    else
                        num_peaks_missed(i_unit) = 0;
                    end
                    
                    num_peaks_amp_assume(end+1,:) = [num_peaks_amp(end,1)+num_peaks_missed(i_unit),num_peaks_amp(end,2)];
                elseif(input_data.amp_list(i_amp)==100) % pass in found value
                    num_peaks_amp_assume(end+1,:) = num_peaks_amp(end,:);
                end
                
            end
            
            
            
        end
        

        
        % get change in latency data for first peak
        data_mask = unit_idx == i_unit & peak_num == 1 & lat_list < 2/1000 & lat_list > 1/1000;
        
        temp_list = lat_list(data_mask) - lat_list(find(data_mask,1,'first'));
        data_mask(find(data_mask,1,'first')) = 0;
        num_add = sum(data_mask);
        
        delta_lat_list(end+1:end+num_add) = temp_list(2:end);
        delta_amp_list(end+1:end+num_add) = amp_list(data_mask);
        delta_unit_idx(end+1:end+num_add) = unit_idx(data_mask);
 
        if(any(temp_list > 0.0005))
            disp('here');
        end
        
    end

    
    
    
    
    % store output_data
    
    output_data.peak_size = peak_size;
    output_data.peak_num = peak_num;
    output_data.lat_list = lat_list;
    output_data.width_list = width_list;
    output_data.std_list = std_list;
    output_data.unit_idx = unit_idx;
    output_data.amp_idx = amp_list;
    
    output_data.num_peaks_amp = num_peaks_amp;
    output_data.num_peaks_amp_assume = num_peaks_amp_assume;
    
    output_data.delta_lat = delta_lat_list;
    output_data.delta_amp = delta_amp_list;
    output_data.delta_unit_idx = delta_unit_idx;
    
    output_data.diam_list = diam_list;
    output_data.cell_list = cell_list;
    output_data.clone_list = clone_list;
    
    output_data.is_model = input_data.is_model;
    
end