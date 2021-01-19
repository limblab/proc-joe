%% load in artifact data
    
    use_duke_board = 1;

    if(use_duke_board)
        load('D:\Lab\Data\stim_ephys_paper\artifact_analysis\sim_study\Duncan_Han_dukeArtifactData.mat');
    else
        load('D:\Lab\Data\stim_ephys_paper\artifact_analysis\sim_study\Duncan_Han_blackrockArtifactData.mat');
    end
    
    load('D:\Lab\Data\stim_ephys_paper\artifact_analysis\sim_study\Neuron_data_lots.mat');

% get stim_on from sync_line_data
    stim_on = cell(size(artifact_data));

    for i_file = 1:numel(sync_line_data)
        stim_on{i_file}=find(diff(sync_line_data{i_file}-mean(sync_line_data{i_file})>3)>.5);
    end


%% remove cases where pulse_widths are not 200us
keep_mask = pulse_width_1==200 & pulse_width_2==200;
amp_1 = amp_1(keep_mask);
amp_2 = amp_2(keep_mask);
artifact_data = artifact_data(keep_mask);
file_list = file_list(keep_mask);
pulse_width_1 = pulse_width_1(keep_mask);
pulse_width_2 = pulse_width_2(keep_mask);
stim_on = stim_on(keep_mask);
sync_line_data = sync_line_data(keep_mask);

%% add neurons to artifact data at various latencies, save NEV file and ground truth data
% use a single stim chan per NEV channel <- will sort multiple channels

    percent_with_spikes = 0.5; % percentage of artifacts with spikes
    base_artifact_idx = 50;
    num_chans_sim = 5;
    num_arts_per_cond = 100;
    
    % pick channels
    unique_chans = unique(stim_chan);
    chan_idx = datasample(1:1:numel(unique_chans),num_chans_sim,'Replace',false);
    
    for i_chan = 1:numel(num_chans_sim)
        % get artifacts related to this channel
        art_data_idx = find(stim_chan == unique_chans(chan_idx(i_chan)));
        amp_list = unique(amp_1(art_data_idx));
        
        for i_art = 1:numel(art_data_idx)
            
        end
    end
    
    
%% full simulation given a method (lol)
% this will randomly give neurons to artifacts at random times, then 






































