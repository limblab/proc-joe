%% load in necessary data
    folderpath = 'C:\Users\jts3256\Desktop\stim_ephys\sim_study\';
    pwd = cd();
    cd(folderpath);
% load in gain as a function of amplitude
    
    load('C:\Users\jts3256\Desktop\stim_ephys\sim_study\duke_board_artifact_gain_ratio.mat');
    waveform_min_peak_idx = 28;
% load in neural data on duke board channel (waveforms) and
% load in artifact data for duke board channel
    amp_list = [5,10,15,20,25,30,40,50,100];
    
    waveforms = []; elec_list = [];
    artifact_data = [];
    amp_num = []; chan_num = [];
    nev_list = dir('nev_data\*.nev');
    output_data_list = dir('nev_data\*outputData*');
    %
    for i_file = 1:numel(nev_list)
        NEV = openNEV([nev_list(i_file).folder,filesep,nev_list(i_file).name],'uV','nosave');
        load([output_data_list(i_file).folder,filesep,output_data_list(i_file).name]);
        
        % map amplitude to amp_list, and remove pw > 200 cases
        wave_keep = ones(numel(outputData.waveforms.parameters),1);
        wave_keep([outputData.waveforms.parameters.pWidth1] > 200 | [outputData.waveforms.parameters.polarity] == 1) = 0;
        wave_keep(any([outputData.waveforms.parameters.amp1]' == amp_list,2)==0) = 0;
        
        amp_mapping = zeros(numel(amp_list),1);
        for i_amp = 1:numel(amp_list)
            if(~isempty(find([outputData.waveforms.parameters.amp1] == amp_list(i_amp) & ...
                wave_keep' == 1)))
                amp_mapping(i_amp) = find([outputData.waveforms.parameters.amp1] == amp_list(i_amp) & ...
                    wave_keep' == 1);
            else
                wave_keep(i_amp) = 0;
            end
        end
        
        % find stims without waveforms within 5 ms
        wave_times = outputData.nevData.ts;
        stim_times = outputData.stimInfo.stimOn/30000; % 30kHz recording
        wave_times = repmat(wave_times',numel(stim_times),1);
        stim_times = repmat(stim_times,1,size(wave_times,2));
        
        time_to_stim = wave_times-stim_times;
        time_to_stim(time_to_stim < 0) = nan;
        time_to_stim = min(time_to_stim',[],'omitnan');
        
        % either no spike at all, or no spike within 6 ms, or max amp from
        % 1.5-6 ms < 100
        filtered_artifact_data = acausalFilter(squeeze(outputData.artifactData.artifact)')';
        window_look = 100 + ceil(30*[1.5,6]);
        stim_mask = (isnan(time_to_stim) | time_to_stim > 0.006) & wave_keep(outputData.waveforms.waveSent)' == 1 & ...
            max(abs(filtered_artifact_data(:,window_look(1):window_look(2))),[],2)' < 40; 
        
        % store stims
        artifact_data = [artifact_data; squeeze(outputData.artifactData.artifact(stim_mask==1,:))];
        amp_num = [amp_num; amp_mapping(outputData.waveforms.waveSent(stim_mask==1))];
        
        % find waveforms on stimulated channel
        str_idx = strfind(nev_list(i_file).name,'Boxchan');
        underscore_idx = strfind(nev_list(i_file).name,'_');
        chan_stim = str2num(nev_list(i_file).name(str_idx+7:underscore_idx(find(underscore_idx > str_idx,1,'first'))-1));
        
        time_to_stim = wave_times-stim_times;
        time_to_stim(time_to_stim < 0) = nan;
        time_to_stim = min(time_to_stim,[],'omitnan');
        wave_mask = NEV.Data.Spikes.Electrode == chan_stim & NEV.Data.Spikes.Unit ~=0 & NEV.Data.Spikes.Unit ~= 255 & (isnan(time_to_stim) | time_to_stim > 0.01);
        temp_wave = mean(NEV.Data.Spikes.Waveform(:,wave_mask)')';
        if(max(temp_wave) < 75*0.254 && max(temp_wave) > 30*0.254)
            waveforms = [waveforms, temp_wave];
            elec_list = [elec_list, chan_stim];
        end
    end
    
    waveforms=waveforms/0.254; 
    
    cd(pwd);
    
%% generate series of artifacts
    
    num_artifacts_per_chan = 100000;
    num_chans = 10;
    prob_spike = 0.5; 
    use_mean_artifact = 0;
    sample_weights = ones(size(artifact_data,1),1);
    
    data = []; wave_data = [];
    
    for i_amp = 1:numel(amp_list)
        sample_weights(amp_num == i_amp) = numel(amp_num)/sum(amp_num==i_amp);
    end
    
    if(use_mean_artifact == 1)
        mean_artifact_data = zeros(numel(amp_list),size(artifact_data,2));
        for i_amp = 1:numel(amp_list)
            art_mask = amp_num == i_amp;
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
            data_temp = acausalFilter(data_temp')';
            data_temp = reshape(data_temp',1,size(data_temp,1)*size(data_temp,2));
            data = [data;data_temp];
            wave_data = [wave_data;amp_num(art_idx)'];
        end
    end
    
    
    
% build gain matrix -- amplifier gain at each data point for each
    % amplitude
    bin_center_idx = floor((bin_centers+0.453)*30);
    
    placed_spike_idx = {};
    placed_elec = {};
    waveform_idx = datasample(1:1:size(waveforms,2),num_chans,'Replace',false);
    for i_chan = 1:num_chans
        artifact_start{i_chan} = 1+((1:num_artifacts_per_chan)-1)*size(artifact_data,2);
        is_spike = rand(num_artifacts_per_chan,1) < prob_spike;
        is_spike(1) = 0; is_spike(end) = 0; % avoid boundary issues
        
        placed_spike_idx{i_chan} = []; placed_elec{i_chan} = [];
        for i_art = 1:num_artifacts_per_chan
            gain_data = ones(size(artifact_data,2),1);

            gain = interp1(bin_center_idx,metric_gain_cathodic(wave_data(i_chan,i_art),:),0:bin_center_idx(end),'linear','extrap');
            gain(isnan(gain)) = 0;
            gain(gain<0) = 0;
            gain(gain>1)=1;
            gain_data(100:100+bin_center_idx(end)) = gain;

            if(is_spike(i_art) == 1)
                idx_post_artifact = ceil(ceil(rand(1)*5*30)+(0.2+0.453)*30); % only place spikes 0.2-7 ms after artifact offset

                placed_spike_idx{i_chan}(end+1,1) = idx_post_artifact + artifact_start{i_chan}(i_art) + 100;
                placed_elec{i_chan}(end+1,1) = elec_list(waveform_idx(i_chan));
                data(i_chan,placed_spike_idx{i_chan}(end)-27:placed_spike_idx{i_chan}(end)+20) = ...
                    data(i_chan,placed_spike_idx{i_chan}(end)-27:placed_spike_idx{i_chan}(end)+20) + ...
                    (gain_data(idx_post_artifact+100-27:idx_post_artifact+100+20).*waveforms(:,waveform_idx(i_chan)))';
            end

        end
    end
            
% get threshold crossings, store waveforms
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
    filename = ['Han_Duncan_20210120_oneWavePerChannel'];
    mapFilename = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    comments = '';
    writeNEV(nevDataAll, packetWidth, filename, mapFilename, comments )
    
%% sort file
%% load in sorted file, then compare threshold crossings to known crossings

    NEV = openNEV(['C:\Users\jts3256\Desktop\stim_ephys\sim_study\',filename,'-s.nev'],'uV');
    
    
%% compare ts to placed_spike_idx
    
    num_bins_combine = 1; % this value times 0.033 ms for each bin size
    % ts is an idx
    tolerance = 15;
    was_placed_spike = {};
    spikes_found = zeros(numel(amp_list),250); % for each amplitude and idx relative to artifact
    spikes_found_ts_list = cell(numel(amp_list),1);
    spikes_total = zeros(size(spikes_found));
    spikes_total_ts_list = cell(numel(amp_list),1);
    
    for i_chan = 1:num_chans
        ts = double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Unit ~= 0 & NEV.Data.Spikes.Unit ~= 255 & NEV.Data.Spikes.Electrode==i_chan));
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
