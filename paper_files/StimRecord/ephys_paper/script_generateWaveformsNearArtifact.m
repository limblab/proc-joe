%% set folder and channel number
folder_path = 'D:\Lab\Data\stim_ephys_paper\long_trains_array_data\example_waves\'; % put \ at end
chan_num = 70;

output_data_files = dir(strcat(folder_path, '*chan', num2str(chan_num), 'stim*outputData*'));
NEV_files = dir(strcat(folder_path, '*chan', num2str(chan_num), 'stim*.nev'));

if(numel(output_data_files) ~= numel(NEV_files))
    error('wrong number of files, please check');
end
%% load in output data and NEV file, grab waveform label from NEV file, put into output data,
%  combine output data's

for i_file = 1:min(4,numel(output_data_files))
    % load files
    load([output_data_files(i_file).folder, filesep, output_data_files(i_file).name]);
    NEV = openNEV([NEV_files(i_file).folder, filesep, NEV_files(i_file).name], 'nosave','uV');
    
    % put waveform label from NEV into outputData
    outputData.rawData.unit = NEV.Data.Spikes.Unit';
    
    % merge output_data's
    if(i_file == 1)
        output_data_all = outputData;
    else
        % all we need is artifactData, rawData, waveforms
        output_data_all.artifactData.t = [output_data_all.artifactData.t; outputData.artifactData.t];
        output_data_all.artifactData.artifact = [output_data_all.artifactData.artifact; outputData.artifactData.artifact];
        
        output_data_all.rawData.ts = [output_data_all.rawData.ts; outputData.rawData.ts];
        output_data_all.rawData.waveforms = [output_data_all.rawData.waveforms; outputData.rawData.waveforms];
        output_data_all.rawData.elec = [output_data_all.rawData.elec; outputData.rawData.elec];
        output_data_all.rawData.unit = [output_data_all.rawData.unit; outputData.rawData.unit];
        
        output_data_all.waveforms.waveSent = [output_data_all.waveforms.waveSent, outputData.waveforms.waveSent];
        
    end
    
    
end

%%  get waveforms for desired unit away from stimulation and immediately after each pulse
unit = 1;
chan_rec = 70;

near_art_window = [1,3]/1000; % in s
far_art_min = 100/1000; % in s

% extract needed data, only use waveforms from the desired unit
artifact_ts = output_data_all.artifactData.t;
unit_ts = output_data_all.rawData.ts(output_data_all.rawData.unit == unit & output_data_all.rawData.elec == chan_rec);
waveforms = output_data_all.rawData.waveforms(output_data_all.rawData.unit == unit & output_data_all.rawData.elec == chan_rec,:);

% get time difference
time_diff = unit_ts - artifact_ts';

% away from artifact if before all artifacts or at least far_art_min after
% an artifact
time_diff_no_neg = time_diff;
time_diff_no_neg(time_diff_no_neg < 0) = 10000; % mask negative times

far_art_mask = all(time_diff < 0, 2) | min(time_diff_no_neg,[],2) > far_art_min;
% near artifact if within near_art_window after an artifact
near_art_mask = min(time_diff_no_neg,[],2) > near_art_window(1) & min(time_diff_no_neg,[],2) < near_art_window(2);

near_waveforms = waveforms(near_art_mask==1,:);
far_waveforms = waveforms(far_art_mask==1,:);

%% align waveforms based on max peak after filtering

min_look = 40;
max_look = 50;
sig_len = 75; % odd number

near_waveforms_filt = acausalFilter(near_waveforms')';
far_waveforms_filt = acausalFilter(far_waveforms')';

[~,near_max_idx] = max(near_waveforms_filt(:,min_look:max_look),[],2);
near_max_idx = near_max_idx + min_look*ones(size(near_waveforms,1),1);

[~,far_max_idx] = max(far_waveforms_filt(:,min_look:max_look),[],2);
far_max_idx = far_max_idx + min_look*ones(size(far_waveforms,1),1);

near_waveforms_adj = zeros(size(near_waveforms,1),sig_len);
far_waveforms_adj = zeros(size(far_waveforms,1),sig_len);

for i = 1:size(near_waveforms_adj,1)
    near_waveforms_adj(i,:) = near_waveforms(i,near_max_idx(i)- (sig_len-1)/2:near_max_idx(i)+(sig_len-1)/2);
end

for i = 1:size(far_waveforms_adj,1)
    far_waveforms_adj(i,:) = far_waveforms(i,far_max_idx(i)- (sig_len-1)/2:far_max_idx(i)+(sig_len-1)/2);
end
%% plot raw far_waveforms, filtered far_waveforms, filtered near_waveforms
max_plot = 50;
far_waveforms_samp = far_waveforms(datasample(1:1:size(far_waveforms,1),max_plot,'Replace',false),:);
near_waveforms_samp = near_waveforms(datasample(1:1:size(near_waveforms,1),max_plot,'Replace',false),:);

figure('Position',[2035 272 255 608]);
ax1=subplot(3,1,1);
x_data = (1:1:size(far_waveforms_samp,2))/30000*1000; 

% raw far waveforms
plot(x_data,far_waveforms_samp - mean(far_waveforms_samp,2),'color',[0.5,0.5,0.5])
formatForLee(gcf);
% filtered far waveforms
ax2=subplot(4,1,2);
plot(x_data,acausalFilter(far_waveforms_samp'),'color',[0.5,0.5,0.5])
formatForLee(gcf);

% raw near waveforms
subplot(4,1,3)
plot(x_data, near_waveforms_samp','color',[0.5,0.5,0.5]);
ylim([100,300])
formatForLee(gcf);
xlim([0,3]);

% filtered near waveforms
ax3=subplot(4,1,4);
plot(x_data,acausalFilter(near_waveforms_samp'),'color',[0.5,0.5,0.5])

linkaxes([ax1,ax2,ax3],'xy');
ylim([-200,200])
xlim([0,3])

formatForLee(gcf);
















