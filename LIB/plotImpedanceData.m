%% load in impedance data somehow (assumes imp is the variable name)

%% remove internal reference and make units kOhms
imp.impedance = imp.impedance(2:end)/1000;

%% make heatmap of impedances
map_filename = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\right S1 20180919\SN 6251-001804.cmp';
max_imp = 200;
min_imp = 0;

impedance_data = zeros(10,10);
map_data = loadMapFile(map_filename);
color_list = inferno();
figure();
% plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5);
hold on
for c = 1:numel(imp.impedance)
    map_idx = find(map_data.chan == c);
    pos = [map_data.row(map_idx),map_data.col(map_idx)];
    impedance_data(pos(1),pos(2)) = imp.impedance(c);
%     rectangle('Position',[pos(1),pos(2),1,1],'EdgeColor',color_list(color_idx,:),'FaceColor',color_list(color_idx,:));
end
colormap(color_list);
imagesc(impedance_data);
b = colorbar;
colormap(color_list);
set(gca,'visible','off')
axis square

