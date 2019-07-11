type = [];
    % 0 =none, 1 = D, 2 = C
    % ascending order, this is for Han currently
    
    
%%
type = [];
for i = 1:numel(P)
    if(strcmpi(P{i},'P'))
        type(i) = 1;
    elseif(strcmpi(P{i},'C'))
        type(i) = 2;
    else
        type(i) = 0;
    end
end
 %% chan_mapping 
 map_data = loadMapFile('R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp')
mapping_data = [map_data.chan,flip(type)'];

%% get RT based on channel number
RT_data = [];
for i = 1:size(mapping_data,1)
    all_files_idx = find(all_files_data.chan == mapping_data(i,1));
    if(~isempty(all_files_idx))
        RT_data(end+1,:) = [all_files_data.mean_rt(all_files_idx),mapping_data(i,2)];
    end
    
end



%% plot cutaneous (2) and deep (1) histograms
data_cut_idx = RT_data(:,2) == 2;
data_deep_idx = RT_data(:,2) == 1;
data_cut = RT_data(data_cut_idx,1);
data_deep = RT_data(data_deep_idx,1);

bin_edges = 0.1:0.025:0.3;
cut_binned = histcounts(data_cut,bin_edges)
deep_binned = histcounts(data_deep,bin_edges)

f = figure();
f.Name = 'Han_singleElectrodes_typeDistributions';

histogram(data_deep,bin_edges,'FaceColor','r')
hold on
histogram(data_cut,bin_edges,'FaceColor','k')
plot(bin_edges(1:end-1)+mode(diff(bin_edges))/2,deep_binned,'r','linewidth',2)
plot(bin_edges(1:end-1)+mode(diff(bin_edges))/2,cut_binned,'k','linewidth',2)

l = legend('Deep','Cutaneous');
set(l,'box','off');
formatForLee(gcf);
set(gca,'fontsize',16);
xlabel('RT (s)');
ylabel('Count');

