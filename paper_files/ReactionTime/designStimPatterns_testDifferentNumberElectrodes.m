% this script assumes that a pattern exists that contains some number of
% electrodes (32 currently). Then this script will build patterns with less
% electrodes which are sampled from the given set

chan_num_all = pattern.chan_num;
max_freq = 330;
num_elec = 16:8:32;

chan_num = {};
wave_freq = {};

[freq_all_norm] = unique(pattern.biomimetic_freq_norm);
freq_all_norm(freq_all_norm == 0) = [];
freq_all = freq_all_norm*max_freq;

freq_all(freq_all < 16) = 16;

chan_num{3} = chan_num_all;
wave_freq{3} = pattern.biomimetic_freq_norm;

keep_idx_24 = [randperm(16,12), randperm(16,12) + 16];
chan_num{2} = chan_num_all(keep_idx_24);
wave_freq{2} = pattern.biomimetic_freq_norm(keep_idx_24);
keep_idx_16 = [randperm(12,8), randperm(12,8) + 8];
chan_num{1} = chan_num{2}(keep_idx_16);
wave_freq{1} = pattern.biomimetic_freq_norm(keep_idx_16);

%
wave_mappings = {};

for num_elec_idx = 1:3
    for bio = 1:2
        wave_mappings{end+1} = zeros(num_elec(num_elec_idx)/2,3); % chan_num, wave_freq, wave_num

        if(bio == 1) % biomimetic
            chan_ = chan_num{num_elec_idx}(1:num_elec(num_elec_idx)/2);
            freq_ = wave_freq{num_elec_idx}(1:num_elec(num_elec_idx)/2);
            for freq_idx = 1:numel(freq_)
                freq_all_idx = find(freq_(freq_idx) == freq_all_norm);
                wave_mappings{end}(freq_idx,2) = freq_all(freq_all_idx);
                wave_mappings{end}(freq_idx,3) = freq_all_idx;
                wave_mappings{end}(freq_idx,1) = chan_(freq_idx);
            end
            
        else % nonbiomimetic
            wave_mappings{end}(:,1) = chan_num{num_elec_idx}(randperm(num_elec(num_elec_idx),num_elec(num_elec_idx)/2))';
            wave_mappings{end}(:,2:3) = wave_mappings{end-1}(:,2:3);
        end
    end
end