idx_repeat = {};
chan_used = [];

for i = 1:numel(all_files_data.chan)
    if(isempty(chan_used) || isempty(find(chan_used == all_files_data.chan(i))))
        chan_idxs = find(all_files_data.chan(i) == all_files_data.chan);
        if(numel(chan_idxs) > 1)
            chan_used(end+1) = all_files_data.chan(i);
            idx_repeat{end+1} = chan_idxs;
        end
    end
end


%% plot groups of repeats

figure();
for i = 1:numel(idx_repeat)
    errorbar(i+linspace(-0.2,0.2,numel(idx_repeat{i})),all_files_data.mean_rt(idx_repeat{i}),...
        all_files_data.std_rt(idx_repeat{i}),'.','markersize',20)
    hold on
end

%% get idea for days between each group
all_files_data.filename(idx_repeat{12})

