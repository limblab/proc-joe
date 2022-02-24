% this scrip expects a folder containing .ns3 (emg data) and .ns5 (stim sync) files
% this goes through a corresponding file, and grabs EMG data around stim.
% It will also let you plot the whole EMG data to look at magnitude during
% movements, etc.


folderpath = 'D:\Lab\Data\StimPDs\Pop\20201023\emg_stim\';

pwd = cd;
cd(folderpath);
NS3_files = dir('*.ns3');
NS5_files = dir('*.ns5');

for i_file = 1:8%:numel(NS3_files)

    % load file
    NS3 = openNSx([folderpath,NS3_files(i_file).name],'uV');
    NS5 = openNSx([folderpath,NS5_files(i_file).name]);


    % get stim on times
    sync_idx = 1;
    num_pulses_per_train = 1;

    NS5_t = ((1:1:NS5.MetaTags.DataPoints)-1)/NS5.MetaTags.SamplingFreq + NS5.MetaTags.Timestamp;
    NS3_t = ((1:1:NS3.MetaTags.DataPoints)-1)/NS3.MetaTags.SamplingFreq + NS3.MetaTags.Timestamp;
    srate = NS3.MetaTags.SamplingFreq;
% 
    stim_on = NS5_t(find(diff(NS5.Data(sync_idx,:)-mean(NS5.Data(sync_idx,:))>3)>.5));
    stim_on = stim_on(1:num_pulses_per_train:end);


    % filter EMG signal
    % high pass at 10 Hz, rectify, low pass at 20 Hz
    [blow,alow] = butter(4,2*20/srate); % butter constructs off 1/2 the sampling frequency!!! 
    [bhigh,ahigh] = butter(4,2*50/srate,'high');
   
    
    NS3.Data = double(NS3.Data);
    for i_emg = 1:size(NS3.Data,1)
        NS3.Data(i_emg,:) = filtfilt(blow,alow,abs(filtfilt(bhigh,ahigh,NS3.Data(i_emg,:))));
    end

    emg_names = deblank({NS3.ElectrodesInfo.Label});
    % get EMG signal around each stim for each channel
    window = [-0.3,0.3]; % s
    window_idx = window*srate;

    EMG_data = nan(numel(stim_on),size(NS3.Data,1),diff(window_idx)+1);

    for i_stim = 1:numel(stim_on)
        stim_on_idx = find(NS3_t > stim_on(i_stim),1,'first');

        EMG_data(i_stim,:,:) = NS3.Data(:,stim_on_idx+window_idx(1):stim_on_idx+window_idx(2));
    end

    mean_EMG_data = squeeze(mean(EMG_data,1));
    baseline_EMG_data = mean(mean_EMG_data(:,1:500),2);


    mean_EMG_data = mean_EMG_data;
    f=figure('Position',[2157,108,560,861]);
    num_emg = 12;
    x_data = (window_idx(1):window_idx(2))/srate;
    for i_emg = 1:num_emg
        subplot(num_emg/2,2,i_emg)
        plot(x_data,mean_EMG_data(i_emg,:)','k','linewidth',2);
    
        ylabel(emg_names{i_emg}(5:end));
        set(gca,'fontsize',12)
        formatForLee(gcf)
        if(i_emg == num_emg)
            xlabel('Time after stim on (s)');
        end
        xlim([x_data(1),x_data(end)])
    end
    f.Name = [NS3_files(i_file).name(1:end-8),'_stimEMG'];
    saveFiguresLIB(f,folderpath,f.Name);
end