%% load in necessary data
    folderpath = 'C:\Users\jts3256\Desktop\stim_ephys\sim_study\';
    pwd = cd();
    cd(folderpath);
% load in waveforms recorded with duke board (same waveforms used for that
% simulation)
    load([folderpath,'duke_recorded_waveforms.mat']);
    
    waveform_min_peak_idx = 28;

% load in artifact data for blackrock stim channel
    amp_list = [5,10,15,20,25,30,40,50,100];
    
    elec_list = [];
    artifact_data = []; artifact_stim_chan = [];
    artifact_amp_idx = []; chan_num = [];
    ns5_list = dir('blackrock_nev_data\*.ns5');
    %
    for i_file = 1:numel(ns5_list)
        disp(ns5_list(i_file).name);
        NSx = openNSx([ns5_list(i_file).folder,filesep,ns5_list(i_file).name],'uV');
        
        % map amplitude to amp_list, and remove pw > 200 cases. Do this
        % based on filename. For these files, we go cathodic, then anodic,
        % then cathodic, then anodic. Ditch anodic examples
        % get pulse widths
        pw1_idx = strfind(ns5_list(i_file).name,'PW1');
        pw2_idx = strfind(ns5_list(i_file).name,'PW2');
        amp1_idx = strfind(ns5_list(i_file).name,'A1');
        amp2_idx = strfind(ns5_list(i_file).name,'A2');
        stim_idx = strfind(ns5_list(i_file).name,'stim');
        chan_idx = strfind(ns5_list(i_file).name,'chan');
        underscore_idx = strfind(ns5_list(i_file).name,'_');
        
        pulse_width_1 = str2num(ns5_list(i_file).name(pw1_idx+4:underscore_idx(find(underscore_idx > pw1_idx,1,'first'))-1));
        pulse_width_2 = str2num(ns5_list(i_file).name(pw2_idx+4:underscore_idx(find(underscore_idx > pw2_idx,1,'first'))-1));
        amp_1 = str2num(ns5_list(i_file).name(amp1_idx+3:underscore_idx(find(underscore_idx > amp1_idx,1,'first'))-1));
        amp_2 = str2num(ns5_list(i_file).name(amp2_idx+3:underscore_idx(find(underscore_idx > amp2_idx,1,'first'))-1));
        stim_chan = str2num(ns5_list(i_file).name(chan_idx+4:stim_idx-1));
        
        % store artifact data into larger list
        if(pulse_width_1 == 200 && pulse_width_2 == 200 && ~isempty(find(amp_list==amp_1)))
            % get stim on times, then artifact data
            NSx_data = double(NSx.Data);
            stim_on=find(diff(NSx_data(end,:)-mean(NSx_data(end,:))>3)>.5);
            
            temp_art_data = [];
            for i_stim = 1:2:numel(stim_on)
                temp_art_data = [temp_art_data; NSx_data(stim_chan,stim_on(i_stim)-100:stim_on(i_stim)+390-1)];
            end
            
            % keep track of stim channel and amplitude
            num_arts = size(temp_art_data,1);
            artifact_data = [artifact_data; temp_art_data];
            artifact_amp_idx = [artifact_amp_idx; find(amp_1==amp_list)+zeros(num_arts,1)];
            artifact_stim_chan = [artifact_stim_chan; stim_chan+zeros(num_arts,1)];
        end
        
    end
    cd(pwd);
    
%% generate series of artifacts
    
    num_artifacts_per_chan = 100000;
    num_chans = 10;
    prob_spike = 0.5; 
    use_mean_artifact = 0;
    sample_weights = ones(size(artifact_data,1),1);
    
    data = []; wave_data = []; is_saturated = [];
    
    for i_amp = 1:numel(amp_list)
        sample_weights(artifact_amp_idx == i_amp) = numel(artifact_amp_idx)/sum(artifact_amp_idx==i_amp);
    end
    
    if(use_mean_artifact == 1)
        mean_artifact_data = zeros(numel(amp_list),size(artifact_data,2));
        for i_amp = 1:numel(amp_list)
            art_mask = artifact_amp_idx == i_amp;
            mean_artifact_data(i_amp,:) = mean(artifact_data(art_mask,:));
        end
        
        art_idx = datasample(1:1:size(mean_artifact_data,1),num_artifacts_per_chan,'Replace',true);
        data = mean_artifact_data(art_idx,:);
        data = acausalFilter(data')';
        data = reshape(data',1,size(data,1)*size(data,2));

        wave_data = art_idx;
    else
        for i_chan = 1:num_chans
            art_idx = datasample(1:1:size(artifact_data,1),num_artifacts_per_chan,'Replace',true,'weights',sample_weights);

            data_temp = artifact_data(art_idx,:);
            is_saturated_temp = abs(data_temp) > 7950;
            data_temp = acausalFilter(data_temp')';
            data_temp = reshape(data_temp',1,size(data_temp,1)*size(data_temp,2));
            is_saturated_temp = reshape(is_saturated_temp',1,size(is_saturated_temp,1)*size(is_saturated_temp,2));
            data = [data;data_temp];
            is_saturated = [is_saturated;is_saturated_temp];
            wave_data = [wave_data;artifact_amp_idx(art_idx)'];
        end
    end
    
% add spike data to data, only if data is not saturated   
    placed_spike_idx = {};
    waveform_idx = datasample(1:1:size(waveforms,2),num_chans,'Replace',false);
    for i_chan = 1:num_chans
        artifact_start{i_chan} = 1+((1:num_artifacts_per_chan)-1)*size(artifact_data,2);
        is_spike = rand(num_artifacts_per_chan,1) < prob_spike;
        is_spike(1) = 0; is_spike(end) = 0; % avoid boundary issues
        
        placed_spike_idx{i_chan} = []; 
        for i_art = 1:num_artifacts_per_chan
            if(is_spike(i_art) == 1)
                idx_post_artifact = ceil(ceil(rand(1)*7.5*30)+(0.2+0.453)*30); % only place spikes 0.2-7 ms after artifact offset

                placed_spike_idx{i_chan}(end+1,1) = idx_post_artifact + artifact_start{i_chan}(i_art) + 100;
                data(i_chan,placed_spike_idx{i_chan}(end)-27:placed_spike_idx{i_chan}(end)+20) = ...
                    data(i_chan,placed_spike_idx{i_chan}(end)-27:placed_spike_idx{i_chan}(end)+20) + ...
                    (waveforms(:,waveform_idx(i_chan)).*(1-is_saturated(i_chan,placed_spike_idx{i_chan}(end)-27:placed_spike_idx{i_chan}(end)+20)'))';
            end

        end
    end
            
%% get threshold crossings, store waveforms
    nevDataAll.ts = []; nevDataAll.waveforms = []; nevDataAll.elec = [];
    for i_chan = 1:num_chans
        threshold = 35;
        data_temp = data(i_chan,:);
        threshold_crossings = find(data_temp > threshold & abs(data_temp) < 1000 & 1:1:numel(data_temp) > 48 & 1:1:numel(data_temp) < numel(data_temp)-48);

        % remove chains -- find largest spot
        idx = 2;
        chain = [1];
        crossings_keep = [];
        while idx <= numel(threshold_crossings)
            if(threshold_crossings(idx) == threshold_crossings(idx-1)+1) % store in chain
                chain = [chain;idx];
            elseif(~isempty(chain)) % broke a chain, store minidx, update idx, empty chain
                [~,max_idx] = max(data_temp(threshold_crossings(chain)));
                if(isempty(crossings_keep))
                    crossings_keep = [threshold_crossings(max_idx+chain(1)-1)];
                else
                    crossings_keep = [crossings_keep;threshold_crossings(max_idx+chain(1)-1)];
                end
                chain = [idx];
            end
            idx = idx+1;
        end
        if(numel(threshold_crossings) > 0)
            threshold_crossings = [crossings_keep;threshold_crossings(end)];
        end

        % go through and weed out ones that are too close to each other
        % prioritize backwards in time
        crossings_mask = ones(numel(threshold_crossings),1);
        for cross = numel(threshold_crossings):-1:2
            if(crossings_mask(cross) == 1) % check time beforehand to see if one is too close
                cross_check = cross-1;
                while cross_check >= 1 && threshold_crossings(cross_check) >= threshold_crossings(cross) - max(30,max(28,20))
                    crossings_mask(cross_check) = 0;
                    cross_check = cross_check-1;
                end
            end
        end  
        threshold_crossings = threshold_crossings(crossings_mask(:) == 1);


    % get waveforms for threshold crossings and store
        spike_waves = zeros(48,numel(threshold_crossings));
        timestamp = zeros(1,numel(threshold_crossings));
        elec_rec = i_chan*ones(1,numel(threshold_crossings)); %ceil(rand(1,numel(threshold_crossings)));
        for i_cross = 1:numel(threshold_crossings)

            spike_waves(:,i_cross) = data_temp(threshold_crossings(i_cross)-22:threshold_crossings(i_cross)+25);
            timestamp(i_cross) = threshold_crossings(i_cross)/30000; % in seconds
        end

        nevDataAll.ts = [nevDataAll.ts;timestamp'];
        nevDataAll.waveforms = [nevDataAll.waveforms;spike_waves'];
        nevDataAll.elec = [nevDataAll.elec;elec_rec'];

    end
    packetWidth = 104;
    cd(folderpath);
    filename = ['Han_Duncan_20210120_oneWavePerChannel_blackrock'];
    mapFilename = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    comments = '';
    writeNEV(nevDataAll, packetWidth, filename, mapFilename, comments )
    
%% sort file
%% load in sorted file, then compare threshold crossings to known crossings

    NSx = openNEV(['C:\Users\jts3256\Desktop\stim_ephys\sim_study\',filename,'-s.nev'],'uV');
    
    
%% compare ts to placed_spike_idx
    
    num_bins_combine = 1; % this value times 0.033 ms for each bin size
    % ts is an idx
    tolerance = 10;
    was_placed_spike = {};
    spikes_found = zeros(numel(amp_list),250); % for each amplitude and idx relative to artifact
    spikes_found_ts_list = cell(numel(amp_list),1);
    spikes_total = zeros(size(spikes_found));
    spikes_total_ts_list = cell(numel(amp_list),1);
    
    for i_chan = 1:num_chans
        ts = double(NSx.Data.Spikes.TimeStamp(NSx.Data.Spikes.Unit ~= 0 & NSx.Data.Spikes.Unit ~= 255 & NSx.Data.Spikes.Electrode==i_chan));
        was_placed_spike{i_chan} = zeros(size(ts));

        for i_spike = 1:numel(ts)
            was_placed_spike{i_chan}(i_spike) = any(abs(ts(i_spike)-placed_spike_idx{i_chan}) <= tolerance);
        end

        % get found spike as a percentage of placed spikes close to artifacts
        % of each stimulation amplitude
        for i_amp = 1:1:numel(amp_list)

            stim_idx = artifact_start{i_chan}(find(wave_data(i_chan,:) == i_amp)) + 100; % 100 data points before stimulation in each chunk

            for i_spike = 1:numel(placed_spike_idx{i_chan})
                % get spikes total and found relative to stim onset. Shift
                % to relative to stim offset when plotting.
                if(any(placed_spike_idx{i_chan}(i_spike) - stim_idx < 8*30 & placed_spike_idx{i_chan}(i_spike) - stim_idx > 0)) % it counts as close to an artifact
                    idx = ceil(min(abs(placed_spike_idx{i_chan}(i_spike)-stim_idx))/num_bins_combine);
                    spikes_total(i_amp,idx) = spikes_total(i_amp,idx) + 1;
                    
                    found_spike = any(abs(placed_spike_idx{i_chan}(i_spike)-ts) <= tolerance);
                    spikes_found(i_amp,idx) = spikes_found(i_amp,idx) + found_spike;
                    spikes_total_ts_list{i_amp}(end+1) = idx;
                    if(found_spike>0)
                        spikes_found_ts_list{i_amp}(end+1) = idx;
                    end
                    
                end
            end
        end
    end

%% plot using spikes_total_ts_list (preferred method)    
    figure();
    amps_plot = [2,5,8,9]; % idk in amp_list
    colors = inferno(numel(amps_plot)+1); % remove most yellow color
    wave_length = 0.453; % ms
    ts_post_offset = 1; % will subtract wave_length from spike times
    bin_size = 0.5; % ms
    
    bin_edges = (0.25:bin_size:7)+wave_length;
    prop_recovered = nan(numel(spikes_found_ts_list),numel(bin_edges)-1);
   
    if(ts_post_offset)
        bin_edges = bin_edges-wave_length;
    end
    
    for i_amp = 1:numel(spikes_found_ts_list)
        % convert to time post stim offset (or onset), then bin and get
        % counts
        temp_found_ts = spikes_found_ts_list{i_amp}/30;
        temp_total_ts = spikes_total_ts_list{i_amp}/30;
        if(ts_post_offset==1)
            temp_found_ts = temp_found_ts-wave_length;
            temp_total_ts = temp_total_ts-wave_length;
        end
        temp_found_counts = histcounts(temp_found_ts,bin_edges);
        temp_total_counts = histcounts(temp_total_ts,bin_edges);
        prop_recovered(i_amp,:) = temp_found_counts./temp_total_counts;
    end
    % normalize prop_recovered by asymptotic value (max in any bin over 4.5ms)
    idx = find(bin_edges > 4.5,1,'first');
    max_recov = max(max(prop_recovered(:,idx:end),[],'omitnan'),[],'omitnan');
    prop_recovered = prop_recovered./max_recov;
    prop_recovered(prop_recovered > 1) = 1;
    
    color_counter = 1;
    for i_amp = amps_plot
        % plot % recovered against idx 
        plot(bin_edges(1:end-1)+mode(diff(bin_edges))/2,prop_recovered(i_amp,:),...
            'color',colors(color_counter,:),'linewidth',1.5)
        hold on
        
        color_counter = color_counter+1;
        xlabel('Time post stim (ms)') 
        ylabel('Proportion of spikes recovered')
    end
    l=legend('10\muA','25\muA','50\muA','100\muA');
    set(l,'box','off','fontsize',14)
    set(gca,'fontsize',14)
    formatForLee(gcf)    



    
    
%% plot using spikes_total (deprecated method)
    figure();
    colors = inferno(7);
    num_bins_combine = 4;
    spikes_found_combine = zeros(size(spikes_found,1),ceil(size(spikes_found,2)/num_bins_combine));
    spikes_total_combine = zeros(size(spikes_found_combine));
    
    for bin_idx = 1:num_bins_combine:size(spikes_found,2)
        spikes_found_combine(:,floor(bin_idx/num_bins_combine)+1) = sum(spikes_found(:,bin_idx:min(bin_idx+num_bins_combine-1,size(spikes_found,2))),2);
        spikes_total_combine(:,floor(bin_idx/num_bins_combine)+1) = sum(spikes_total(:,bin_idx:min(bin_idx+num_bins_combine-1,size(spikes_total,2))),2);
    end
    
    color_counter = 1;
    for i_amp = [2,4,7,9]
        % plot % recovered against idx 
        plot_mask = ~isnan(spikes_total_combine(i_amp,:));
        plot(0.033*num_bins_combine*(0:1:size(spikes_total_combine,2)-1) + 0.033*num_bins_combine/2 - 0.453,...
            spikes_found_combine(i_amp,plot_mask)./spikes_total_combine(i_amp,plot_mask),...
            'color',colors(color_counter,:),'linewidth',1.5)
        hold on
        
        color_counter = color_counter+1;
        xlabel('Time post stim (ms)') 
        ylabel('Proportion of spikes recovered')
    end
%     l=legend('5\muA','10\muA','20\muA','30\muA','50\muA','100\muA');
%     set(l,'box','off','fontsize',14)
    set(gca,'fontsize',14)
    formatForLee(gcf)
