%% do a bar plot
counter = 1;
colors = {'r','b',[204,102,0]/255};
num_elecs = [3,6,12];
total_charge = [240,360,480];
offset = [-5,0,5];
bar_data = [];
figure();
hold on
for num_elec_idx = 1:numel(num_elecs)
    for total_charge_idx = 1:numel(total_charge)
        mean_rt = mean(data.cueInfo(counter).rt);
        std_rt = std(data.cueInfo(counter).rt);
        bar_data(counter).num_elecs = num_elecs(num_elec_idx);
        bar_data(counter).total_charge = total_charge(total_charge_idx);
        bar_data(counter).y_mean = mean_rt;
        bar_data(counter).y_std = std_rt;
        counter = counter + 1;
    end
end

for counter = 1:numel(bar_data)
    plot(bar_data(counter).total_charge+offset(find(num_elecs==bar_data(counter).num_elecs)),...
        bar_data(counter).y_mean,'.','color',colors{find(num_elecs==bar_data(counter).num_elecs)},'markersize',16)
end

%% 
f = gcf;
f.Name = 'Han_20180807_multielec_rtVsNumChannels';
saveFiguresLIB(f,inputData.folderpath,f.Name);

%% old, dots instead of a bar
counter = 1;
colors = {'r','b',[204,102,0]/255};
num_elecs = [3,6,12];
total_charge = [240,360,480];
offset = [-0.025,0,0.025];
figure();
hold on
for num_elec_idx = 1:numel(num_elecs)
    for total_charge_idx = 1:numel(total_charge)
        mean_rt = mean(data.cueInfo(counter).rt);
        std_rt = std(data.cueInfo(counter).rt);
        plot(log(num_elecs(num_elec_idx))/log(3)+offset(total_charge_idx),mean_rt,'.','color',colors{total_charge_idx},'markersize',16)
        counter = counter + 1;
    end
end
counter = 1;
for num_elec_idx = 1:numel(num_elecs)
    for total_charge_idx = 1:numel(total_charge)
        mean_rt = mean(data.cueInfo(counter).rt);
        std_rt = std(data.cueInfo(counter).rt);
        plot(offset(total_charge_idx)+log([num_elecs(num_elec_idx),num_elecs(num_elec_idx)])/log(3),...
            mean_rt+std_rt*[-1,1],'-','color',colors{total_charge_idx},'linewidth',1.5)
        counter = counter + 1;
    end
end

mean_bump_rt = mean(data.cueInfo(end).rt);
std_bump_rt = std(data.cueInfo(end).rt);
x_data = log(num_elecs(1:2:3)) + [-0.5,0.2];
plot(x_data,[mean_bump_rt,mean_bump_rt],'k','linewidth',1.5)
plot(x_data,[mean_bump_rt+std_bump_rt,mean_bump_rt+std_bump_rt],'k--','linewidth',1.5)
plot(x_data,[mean_bump_rt-std_bump_rt,mean_bump_rt-std_bump_rt],'k--','linewidth',1.5)

xlim([x_data])
ylim([0.1,0.3])
ax = gca;
ax.XTick = log(num_elecs)/log(3);
ax.XTickLabel = num2str(num_elecs');
ax.FontSize = 14;
xlabel('Num stimulation channels');
ylabel('RT (s)');
formatForLee(gcf)
ax.XMinorTick = 'off';
l = legend('240\muA','360','480');
set(l,'box','off')
