%% setup folders
array_data_folder = 'E:\Data\Joseph\long_trains_array_data\amp_freq\Han_adj\';
NS5_folder = 'E:\Data\Joseph\long_trains_array_data\amp_freq\NS5\';
array_data_files = dir([array_data_folder,'*arrayData.mat']);

%% insert stim on times

% for each array data in array_data_folder:
for array_file_idx = 1:numel(array_data_files)
    % load arrayData
    load([array_data_folder,array_data_files(array_file_idx).name]);
    
    % find corresponding NS5 files based on stim chan
    stim_chan = num2str(arrayData{1}.CHAN_LIST{1});
    NS5_data_files = dir([NS5_folder,'*chan',stim_chan,'stim*','spikesExtracted.ns5*']);
    % these are in numerical (chronological) order
    
    % for each file, store when each pulse happened per condition.
    
    pulse_time = cell(numel(arrayData{1}.binCounts),1); % time of pulse for each condition
    wave_counter = 1;
    for NS5_idx = 1:numel(NS5_data_files)
        % open NS5 file
        NS5 = openNSx([NS5_folder,NS5_data_files(NS5_idx).name],'uV');
        sync_data = NS5.Data(2,:);
        time_data = (1:1:numel(sync_data))/30000;
        % find all stim on's
        stim_on = time_data(find(diff(sync_data-mean(sync_data)>3)>.5));
        
        % find start of trains
        idx_start = find(diff(stim_on) > 1);
        idx_start = [1,idx_start + 1];
        train_on = stim_on(idx_start);
        
        % for each train, store pulse times into pulse_time
        for t = 1:numel(train_on)
            pulse_train = stim_on - train_on(t);
            keep_mask = pulse_train >= 0 & pulse_train < 4.5;
            pulse_train = pulse_train(keep_mask);
            
            pulse_time{arrayData{1}.WAVEFORM_SENT(wave_counter)}{end+1} = pulse_train;
            wave_counter = wave_counter + 1;
        end
    end
    
    % save array data with updated filename
    for unit = 1:numel(arrayData)
        arrayData{unit}.PULSE_TIMES = pulse_time;
    end
    save([array_data_folder,array_data_files(array_file_idx).name(1:end-4),'_pulseTimes'],'arrayData','inputData','optsExtract');
    
end