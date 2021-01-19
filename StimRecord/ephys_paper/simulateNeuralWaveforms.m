%% load in necessary data
    folderpath = 'E:\Data\Joseph\sim_neural_data\';
    pwd = cd();
    cd(folderpath);
% load in gain as a function of amplitude
    
    load('E:\Data\Joseph\sim_neural_data\duke_board_artifact_gain_ratio.mat');
    waveform_min_peak_idx = 28;
    CHAN_WAVEFORM = 95;
% load in neural data on duke board channel (waveforms) and
% load in artifact data for duke board channel
    amp_list = [5,10,15,20,25,30,40,50,100];
    
    waveforms = []; elec_list = [];
    artifact_data = [];
    amp_num = [];
    nev_list = dir('*.nev');
    output_data_list = dir('*outputData*');
    for i_file = 1:numel(nev_list)
        NEV = openNEV([folderpath,nev_list(i_file).name],'uV','nosave');
        load([folderpath,output_data_list(i_file).name]);
        
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
        waveforms = [waveforms, mean(NEV.Data.Spikes.Waveform(:,wave_mask)')'];
        elec_list = [elec_list, chan_stim+zeros(1,sum(wave_mask))];
    end
    
    waveforms=waveforms/0.254; 
    cd(pwd);
    
%% generate series of artifacts
    
    num_artifacts = 400000;
    prob_spike = 0.6; 
    use_mean_artifact = 0;
    sample_weights = ones(size(artifact_data,1),1);
    for i_amp = 1:numel(amp_list)
        sample_weights(amp_num == i_amp) = numel(amp_num)/sum(amp_num==i_amp);
    end
    
    if(use_mean_artifact == 1)
        mean_artifact_data = zeros(numel(amp_list),size(artifact_data,2));
        for i_amp = 1:numel(amp_list)
            art_mask = amp_num == i_amp;
            mean_artifact_data(i_amp,:) = mean(artifact_data(art_mask,:));
        end
        
        art_idx = datasample(1:1:size(mean_artifact_data,1),num_artifacts,'Replace',true);
        data = mean_artifact_data(art_idx,:);
        data = acausalFilter(data')';
        data = reshape(data',1,size(data,1)*size(data,2));

        wave_data = art_idx;
    else
        art_idx = datasample(1:1:size(artifact_data,1),num_artifacts,'Replace',true,'weights',sample_weights);
    
        data = artifact_data(art_idx,:);
        data = acausalFilter(data')';
        data = reshape(data',1,size(data,1)*size(data,2));
        wave_data = amp_num(art_idx);
    end
    
    
    artifact_start = 1+((1:numel(wave_data))-1)*size(artifact_data,2);
    % build gain matrix -- amplifier gain at each data point for each
    % amplitude
    bin_center_idx = floor((bin_centers+0.453)*30);
    is_spike = rand(num_artifacts,1) < prob_spike;
    is_spike(1) = 0; is_spike(end) = 0; % avoid boundary issues
    placed_spike_idx = [];
    placed_elec = [];
    for i_art = 1:numel(art_idx)
        gain_data = ones(size(artifact_data,2),1);

        gain = interp1(bin_center_idx,median_gain_cathodic(wave_data(i_art),:),bin_center_idx(1):bin_center_idx(end));
        
        gain_data(100+bin_center_idx(1)-1:100+bin_center_idx(end)-1) = gain;
        
        if(is_spike(i_art) == 1)
            waveform_idx = datasample(1:1:size(waveforms,2),1);
            idx_post_artifact = ceil(ceil(rand(1)*6*30)+(0.5+0.453)*30); % only place spikes 0.5-6 ms after artifact onset
            
            placed_spike_idx(end+1,1) = idx_post_artifact + artifact_start(i_art) + 100;
            placed_elec(end+1,1) = elec_list(waveform_idx);
            data(placed_spike_idx(end)-27:placed_spike_idx(end)+20) = data(placed_spike_idx(end)-27:placed_spike_idx(end)+20) + ...
                (gain_data(idx_post_artifact+100-27:idx_post_artifact+100+20).*waveforms(:,waveform_idx))';
        end
        
    end
            
%% get threshold crossings, store waveforms

    threshold = 40;
    threshold_crossings = find(data > threshold & abs(data) < 1000 & 1:1:numel(data) > 48 & 1:1:numel(data) < numel(data)-48);
    
    % remove chains -- find largest spot
    idx = 2;
    chain = [1];
    crossings_keep = [];
    while idx <= numel(threshold_crossings)
        if(threshold_crossings(idx) == threshold_crossings(idx-1)+1) % store in chain
            chain = [chain;idx];
        elseif(~isempty(chain)) % broke a chain, store minidx, update idx, empty chain
            [~,max_idx] = max(data(threshold_crossings(chain)));
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
    elec_rec = ones(1,numel(threshold_crossings)); %ceil(rand(1,numel(threshold_crossings)));
    for i_cross = 1:numel(threshold_crossings)
        
        spike_waves(:,i_cross) = data(threshold_crossings(i_cross)-22:threshold_crossings(i_cross)+25);
        timestamp(i_cross) = threshold_crossings(i_cross)/30000; % in seconds
    end

    packetWidth = 104;
    filename = ['test'];
    nevDataAll.ts = timestamp';
    nevDataAll.waveforms = spike_waves';
    nevDataAll.elec = elec_rec';
    mapFilename = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    comments = '';
    writeNEV(nevDataAll, packetWidth, filename, mapFilename, comments )
    
%% sort file
%% load in sorted file, then compare threshold crossings to known crossings

    NEV = openNEV('C:\Users\jts3256\Desktop\GIT\test-01.nev','uV');
    
    
%% compare ts to placed_spike_idx
    ts = double(NEV.Data.Spikes.TimeStamp(NEV.Data.Spikes.Unit ~= 0 & NEV.Data.Spikes.Unit ~= 255));
    num_bins_combine = 1; % this value times 0.033 ms for each bin size
    % ts is an idx
    tolerance = 10;
    was_placed_spike = zeros(size(ts));
    
    for i_spike = 1:numel(ts)
        was_placed_spike(i_spike) = any(abs(ts(i_spike)-placed_spike_idx) <= tolerance);
    end
    
    % get found spike as a percentage of placed spikes close to artifacts
    % of each stimulation amplitude
    spikes_found = zeros(size(median_gain_cathodic,1),220); % for each amplitude and idx relative to artifact
    spikes_total = zeros(size(median_gain_cathodic,1),220);
    
    for i_amp = 1:1:9%:size(median_gain_cathodic,1)
        
        stim_idx = artifact_start(find(wave_data == i_amp)) + 100; % 100 data points before stimulation in each chunk
        
        for i_spike = 1:numel(placed_spike_idx)
            if(any(placed_spike_idx(i_spike) - stim_idx < 7*30 & placed_spike_idx(i_spike) - stim_idx > 0)) % it counts as close to an artifact
                idx = ceil(min(abs(placed_spike_idx(i_spike)-stim_idx))/num_bins_combine);
                spikes_total(i_amp,idx) = spikes_total(i_amp,idx) + 1;
                spikes_found(i_amp,idx) = spikes_found(i_amp,idx) + any(abs(placed_spike_idx(i_spike)-ts) <= tolerance);
            end
        end
   end

    %%
    figure();
    colors = inferno(7);
    num_bins_combine = 15;
    spikes_found_combine = zeros(size(spikes_found,1),ceil(size(spikes_found,2)/num_bins_combine));
    spikes_total_combine = zeros(size(spikes_found_combine));
    
    for bin_idx = 1:num_bins_combine:size(spikes_found,2)
        spikes_found_combine(:,floor(bin_idx/num_bins_combine)+1) = sum(spikes_found(:,bin_idx:min(bin_idx+num_bins_combine-1,size(spikes_found,2))),2);
        spikes_total_combine(:,floor(bin_idx/num_bins_combine)+1) = sum(spikes_total(:,bin_idx:min(bin_idx+num_bins_combine-1,size(spikes_total,2))),2);
    end
    
    color_counter = 1;
    for i_amp = [1,2,4,6,8,9]
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
    l=legend('5\muA','10\muA','20\muA','30\muA','50\muA','100\muA');
    set(l,'box','off','fontsize',14)
    set(gca,'fontsize',14)
    formatForLee(gcf)
%% find time when amp crosses threshold and plot that
    thresh = 0.7*0.75;
    
    idx_cross = zeros(numel(amp_list),1);
    
    figure();
    for i_amp = 1:numel(amp_list)
        plot_mask = ~isnan(spikes_total_combine(i_amp,:));
        idx_cross(i_amp) = find(spikes_found_combine(i_amp,plot_mask)./spikes_total_combine(i_amp,plot_mask) > thresh,1,'first');
    end
    
    plot(amp_list,idx_cross)
    
    %%
    figure();
    x_data = ((1:1:size(artifact_data,2))-101)/30 - 0.453;
    amp = 100;
    amp_idx = find(amp_list == amp);
    
    art_idx = find(amp_num == amp_idx);
    
    plot(x_data,mean(acausalFilter(artifact_data(art_idx,:)')'))
