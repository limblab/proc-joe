plotted_data_temp = [];
best_r2 = 0;
best_sub = 0;
all_r2 = [];
for i = 0
    plotted_data_temp = plotted_data;
    plotted_data_temp(:,1) = plotted_data_temp(:,1) - i;
    r2 = 1 - sum((plotted_data_temp(:,2)-plotted_data_temp(:,1)).^2)/sum((plotted_data_temp(:,2)-mean(plotted_data_temp(:,2))).^2);
    all_r2(end+1) = r2;
    if(r2 > best_r2)
        best_r2 = r2;
        best_sub = i;
    end
end
%1 - sum((plotted_data(:,2)-plotted_data(:,1)).^2)/sum((plotted_data(:,2)-mean(plotted_data(:,2))).^2)