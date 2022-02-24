folderpath = 'D:\Lab\Data\StimEvokedMove\Han_20200910_emg\';
pwd = cd;
cd(folderpath);
filename = dir('*.ns5');

openNSx([folderpath,filename(1).name],'uV');

cd(pwd)

emg_idx = 2;
stim_idx = 3;
num_stims_per_cond = 20;
num_data_pre = 6000;
num_data_post = 12000;
sample_rate = 30000; % Hz

emg_data_all = NS5.Data(emg_idx,:);
sync_data = NS5.Data(stim_idx,:);

stim_on = find(diff(sync_data-mean(sync_data)>3)>.5);
stim_off = find(diff(sync_data-mean(sync_data)<3)>.5);



%% bandpass filter (80-500 Hz) between stimulation pulses. Blank during stim
% rectify
% low pass at 10Hz
data_post_stim_blank = 60;

[bpass_filt_b,bpass_filt_a] = butter(1,[80,500]/(sample_rate/2),'bandpass');
[low_b,low_a] = butter(1,[10]/(sample_rate/2),'low');

% filter before first pulse
emg_data_filt = emg_data_all;
median_emg = median(emg_data_all);
emg_data_filt(1:stim_on(1)-1) = acausalFilter(emg_data_filt(1:stim_on(1)-1)',bpass_filt_b,bpass_filt_a)';
emg_data_filt(1:stim_on(1)-1) = abs(emg_data_filt(1:stim_on(1)-1));
emg_data_filt(1:stim_on(1)-1) = acausalFilter(emg_data_filt(1:stim_on(1)-1)',low_b,low_a)';


% filter between each pulse and blank during each pulse (rectify and low
% pass)
for i_stim = 1:numel(stim_on)-1
    emg_data_filt(stim_on(i_stim):stim_off(i_stim)+data_post_stim_blank-1) = 0;
    emg_data_filt(stim_off(i_stim)+data_post_stim_blank:stim_on(i_stim+1)) = acausalFilter(emg_data_filt(stim_off(i_stim)+data_post_stim_blank:stim_on(i_stim+1))',bpass_filt_b,bpass_filt_a)';
    emg_data_filt(stim_off(i_stim)+data_post_stim_blank:stim_on(i_stim+1)) = abs(emg_data_filt(stim_off(i_stim)+data_post_stim_blank:stim_on(i_stim+1)));
    emg_data_filt(stim_off(i_stim)+data_post_stim_blank:stim_on(i_stim+1)) = acausalFilter(emg_data_filt(stim_off(i_stim)+data_post_stim_blank:stim_on(i_stim+1))',low_b,low_a)';
end
% filter after last pulse
emg_data_filt(stim_on(end):stim_off(end)+data_post_stim_blank-1) = 0;
emg_data_filt(stim_off(end)+data_post_stim_blank:end) = acausalFilter(emg_data_filt(stim_off(end)+data_post_stim_blank:end)',bpass_filt_b,bpass_filt_a)';
emg_data_filt(stim_off(end)+data_post_stim_blank:end) = abs(emg_data_filt(stim_off(end)+data_post_stim_blank:end));
emg_data_filt(stim_off(end)+data_post_stim_blank:end) = acausalFilter(emg_data_filt(stim_off(end)+data_post_stim_blank:end)',low_b,low_a)';

emg_near_stim = zeros(numel(stim_on),num_data_pre+num_data_post+1);

for i_stim = 1:numel(stim_on)
    emg_near_stim(i_stim,:) = emg_data_filt(stim_on(i_stim)-num_data_pre:stim_on(i_stim)+num_data_post);
end

%%
cond_map = [1:numel(stim_on)/num_stims_per_cond];
cond_map = repmat(cond_map,num_stims_per_cond,1);
cond_map = reshape(cond_map,numel(cond_map),1);

conds = unique(cond_map);

for i_cond = 1:numel(conds)
    figure();
    
end



