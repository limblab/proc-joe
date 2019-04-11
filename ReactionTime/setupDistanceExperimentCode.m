% mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
mapData = loadMapFile(mapFileName);



chans_pos = [];

chans_to_use = 1:96;
chans_to_use(all_files_data.no_response) = [];

for i = 1:numel(chans_to_use)
    idx = find(all_files_data.chan == chans_to_use(i));
    [rt_use(i)] = min(all_files_data.mean_rt(idx));
    
    idx = find(mapData.chan == chans_to_use(i));
    chans_pos(i,:) = [mapData.row(idx),mapData.col(idx)];
end

%%
f = figure();
b = colorbar;
image_data = zeros(10,10)+0.3;
for i = 1:size(chans_pos,1)
    image_data(chans_pos(i,1),chans_pos(i,2)) = rt_use(i);
end
colormap(flip(inferno))
imagesc(image_data)
for i = 1:size(chans_pos,1)
    hold on
    text(chans_pos(i,2),chans_pos(i,1),num2str(chans_to_use(i)),'color',[1,1,1]);
end
b = colorbar;

%%
chan_num = 88;
chan_idx = find(chans_to_use == chan_num);

other_chans = chans_to_use;
other_chans_rt = rt_use;
other_chans(chan_idx) = [];
other_chans_rt(chan_idx) = [];

rt_diff = abs(rt_use(chan_idx) - other_chans_rt);
[rt_diff_sort,sort_idx] = sort(rt_diff);

rt_use(chan_idx)
[other_chans(sort_idx(1:8))',other_chans_rt(sort_idx(1:8))']

%%
list1 = [68,66,33]; list2 = [68,77,47];
for i = 1:numel(list1)
    idx1(i) = find(list1(i) == chans_to_use);
    idx2(i) = find(list2(i) == chans_to_use);
end

[rt_use(idx1)',rt_use(idx2)']


%%
figure
hold on
x_vals = [1,1,2,2,3,3,4,4,5,5,6,6,7,7] + repmat([-0.25,0.25],1,7)
for i = 1:8%numel(electrodeList)
    for j = 1:numel(electrodeList{i})
        plot(x_vals(i),all_files_data.mean_rt(find(all_files_data.chan == electrodeList{i}(j))),...
            'd','markersize',5,'color',colorList{i},'markerfacecolor',colorList{i});
        hold on
    end
    
end