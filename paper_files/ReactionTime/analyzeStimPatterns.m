%% analyses stim patterns provided reaction time data

% reaction time data comes in a struct from analyze_reactionTime script
% wave_mappings is the corresponding mat file used to build the stim
% pattern

% this script does random things to try and explain the reaction time
% results. woo

%% some setup
cue_info_idx = [];
mean_stim_rt = [];
std_stim_rt = [];
for stim_code = 1:numel(wave_mappings)
    cue_info_idx(stim_code) = find([data.cueInfo.stimCode] == stim_code-1 & [data.cueInfo.bumpMag] == 0);
    mean_stim_rt(stim_code) = mean(data.cueInfo(cue_info_idx(stim_code)).rt);
    std_stim_rt(stim_code) = std(data.cueInfo(cue_info_idx(stim_code)).rt);
end

mean_bump_rt = mean(data.cueInfo(end).rt);
std_bump_rt = std(data.cueInfo(end).rt);
biomimetic_idx = 1:2:numel(mean_stim_rt);
nonbiomimetic_idx = 2:2:numel(mean_stim_rt);
%% compute a 'charge-delivered' for each waveform, plot RT against charge delivered
% since we have the same amplitude, and pulse width, this is really
% just a summation of the frequency for each pattern

charge_delivered = zeros(numel(wave_mappings),1);
for i = 1:numel(wave_mappings)
    charge_delivered(i) = sum(wave_mappings{i}(:,2));
end

figure()
plot(charge_delivered(biomimetic_idx),mean_stim_rt(biomimetic_idx),'k.','markersize',12)
hold on
plot(charge_delivered(nonbiomimetic_idx),mean_stim_rt(nonbiomimetic_idx),'r.','markersize',12)
xlabel('Pseudo charge delivered')
ylabel('RT (s)')
formatForLee(gcf)
l=legend('biomimetic','nonbiomimetic');
set(l,'box','off')
set(gca,'fontsize',14);


%% check to see how 'biomimetic' each mapping is
% plot a circle plot, vector for each channel, dir of vector is PD of
% channel, length of vector is mag of stimulation


for i = 1:numel(wave_mappings)
    figure();
    
    for elec = 1:size(wave_mappings{i},1)
        pd_idx = find(wave_mappings{i}(elec,1) == input_data.chan_list);
        pd = input_data.pd_elec(pd_idx);
        freq_norm = wave_mappings{i}(elec,2);
        polarplot([0,pd],[0,freq_norm],'k-','linewidth',1.5);
        hold on
    end
    
%     polarplot([0,pattern{ceil(i/2)}.dir],[0,1],'r-','linewidth',1.5);
    polarplot([0,input_data.dir_model_all(i)],[0,1],'r-','linewidth',1.5);
    f=gcf;
%     f.Name = ['Han_20180725_rt_polar_',num2str(i)];
%     saveFiguresLIB(f,inputData.folderpath,f.Name);
end

%% draw a heatmap for each mapping
map_filename = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
map_data = loadMapFile(map_filename);
colors = [(0:1:63)'/63 zeros(64,1) zeros(64,1)];

for i = 1:numel(wave_mappings)
    figure();
    plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5)
    hold on
    for elec = 1:size(wave_mappings{i},1)
        map_idx = find(wave_mappings{i}(elec,1) == map_data.chan);
        freq_norm = wave_mappings{i}(elec,2);
        color_idx = max(1,min(size(colors,1),ceil(freq_norm*size(colors,1))));
        
        rectangle('Position',[map_data.row(map_idx),map_data.col(map_idx),1,1],'FaceColor',colors(color_idx,:), 'EdgeColor',colors(color_idx,:));
        hold on
    end
    
    f=gcf;
    set(gca,'visible','off')
    f.Name = ['Han_20181030_rt_heatmap_',num2str(i)];
%     saveFiguresLIB(f,inputData.folderpath,f.Name);
end

%% make wave_mappings charge conserved
total_charge = sum(wave_mappings{1}(:,2));
charge_conserved_wave_mappings = wave_mappings;
for w_idx = 1:numel(wave_mappings)
    charge_conserved_wave_mappings{w_idx}(:,2) = total_charge*charge_conserved_wave_mappings{w_idx}(:,2)/sum(charge_conserved_wave_mappings{w_idx}(:,2));
end


%% dissimilarity metric



