% downsample analog tables to DLC time
cds = cds_list{2};
sync_ana_idx = 2;
sync_name = 'expsync';
dlc_table_idx = 3;

exp_sync = [cds.analog{sync_ana_idx}.t,cds.analog{sync_ana_idx}.(sync_name)];
dlc_tab = cds.analog{dlc_table_idx};

exp_sync = exp_sync(1:100:end,:);
is_exp_sync = exp_sync(:,2) > 2500;

%% label DLC table as in experiment or not (based on button).
is_exp_dlc = zeros(size(dlc_tab,1),1);

for i = 1:numel(is_exp_dlc)
    % find closest sync point, grab label
    t_dlc = dlc_tab.t(i);
    t_dist = exp_sync(:,1) - t_dlc;
    [~,min_idx] = min(abs(t_dist));
    is_exp_dlc(i) = is_exp_sync(min_idx);
end

%% remove entries if elbow1 and hand2 are missing
dlc_tab_names = dlc_tab.Properties.VariableNames;

dlc_idx = [find((strcmpi(dlc_tab_names,['hand2','_x']))),...
                            find((strcmpi(dlc_tab_names,['hand2','_y']))),...
                            find((strcmpi(dlc_tab_names,['hand2','_z']))),...
                            find((strcmpi(dlc_tab_names,['elbow1','_x']))),...
                            find((strcmpi(dlc_tab_names,['elbow1','_y']))),...
                            find((strcmpi(dlc_tab_names,['elbow1','_y'])))];
 
dlc_mat = table2array(dlc_tab);
bad_points = any(isnan(dlc_mat(:,dlc_idx)),2);


keep_mask = ~bad_points & is_exp_dlc;
% compute correlations


%%
    corr_markernames = {'elbow1','hand2'};
    dlc_idx = [];
    for i_marker = 1:numel(corr_markernames)
        dlc_idx = [dlc_idx,find((strcmpi(dlc_tab_names,[corr_markernames{i_marker},'_x']))),...
            find((strcmpi(dlc_tab_names,[corr_markernames{i_marker},'_y']))),...
            ];
    end

    % get velocity for each marker
    dlc_vel = dlc_mat;
    for i = 1:size(dlc_mat,2)
        dlc_vel(:,i) = gradient(dlc_mat(:,i),0.05);
    end
    
    dlc_data_vel = dlc_vel(keep_mask==1,dlc_idx);
    % remove nan's
    dlc_data_vel = dlc_data_vel(~any(isnan(dlc_data_vel),2),:);
    %Calculate speed
    dlc_data_speed = [sqrt(sum(dlc_data_vel(:,1:2).^2,2)),sqrt(sum(dlc_data_vel(:,3:4).^2,2))];        

    corr_vel = corr(dlc_data_vel);
    corr_speed = corr(dlc_data_speed);