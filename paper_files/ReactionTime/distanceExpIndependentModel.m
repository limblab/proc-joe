%% make relevant data
% need independent mean and std, then combined mean and std for all
% pairs/triplets

% so load in data_all from distance experiment (made w/ distanceExpStats.m)

monkey_name = 'Duncan';

%% get mean RT for each multistim pair and store

    independent_model.multi_mean_rt = [];
    independent_model.multi_std_rt = [];
    independent_model.multi_std_err_rt = [];

    independent_model.single_mean_rt = {};
    independent_model.single_std_rt = {};
    independent_model.single_num_trials = [];
    independent_model.is_adjacent = [];
    independent_model.group = [];
    independent_model.num_chans = [];
    
for group = unique(data_all.group) % for each pair of stim groups
    mean_near = mean(data_all.rt(data_all.group==group & data_all.dist==1));
    std_near = std(data_all.rt(data_all.group==group & data_all.dist==1));
    std_err_near = std_near./sqrt(numel(data_all.rt(data_all.group==group & data_all.dist == 1)));
    mean_far = mean(data_all.rt(data_all.group==group & data_all.dist==0));
    std_far = std(data_all.rt(data_all.group==group & data_all.dist==0));
    std_err_far = std_near./sqrt(numel(data_all.rt(data_all.group==group & data_all.dist == 0)));
    
    independent_model.multi_mean_rt(end+1) = mean_near;
    independent_model.multi_std_rt(end+1) = std_near;
    independent_model.multi_std_err_rt(end+1) = std_err_near;
    independent_model.is_adjacent(end+1) = 1;
    independent_model.group(end+1) = group;
    
    independent_model.multi_mean_rt(end+1) = mean_far;
    independent_model.multi_std_rt(end+1) = std_far;
    independent_model.multi_std_err_rt(end+1) = std_err_far;
    independent_model.is_adjacent(end+1) = 0;
    independent_model.group(end+1) = group;
    
    % get single_mean_rt and std_rt for each channel
    
    for i = [1,0] % near, then far
        chan_list = data_all.chans{find(data_all.group == group & data_all.dist == i,1,'first')};
        single_mean_rts = [];
        single_std_rts = [];
        single_num_trials = [];
        for chan = chan_list
            all_files_idx = find(all_files_data.chan == chan);
            single_mean_rts(end+1) = all_files_data.mean_rt(all_files_idx);
            single_std_rts(end+1) = all_files_data.std_rt(all_files_idx);
            single_num_trials(end+1) = all_files_data.num_trials(all_files_idx);
        end
        independent_model.single_mean_rt{end+1} = single_mean_rts;
        independent_model.single_std_rt{end+1} = single_std_rts;
        independent_model.single_num_trials{end+1} = single_num_trials;
        independent_model.num_chans(end+1) = numel(independent_model.single_mean_rt{end});

    end
    
end


%% for each electrode group, get mean and std rt using the race model -- bootstrap std. err of mean
num_trials = 25;
independent_model.race_model_mean = [];
independent_model.race_model_std = [];

for i = 1:numel(independent_model.multi_mean_rt)
    trial_rt = zeros(num_trials,numel(independent_model.single_mean_rt{i}));
    for elec = 1:numel(independent_model.single_mean_rt{i})
        trial_rt(:,elec) = normrnd(independent_model.single_mean_rt{i}(elec),independent_model.single_std_rt{i}(elec),num_trials,1);
    end
    
    min_trial_rt = min(trial_rt,[],2);
    independent_model.race_model_mean(i) = mean (min_trial_rt);
    independent_model.race_model_std(i) = std(min_trial_rt);
end

% compare race model mean with actual model
f=figure();
f.Name = [monkey_name,'_raceModelPrediction'];

% legend placeholders
plot(-100,-100,'k.','markersize',20)
hold on
plot(-100,-100,'ko','markersize',10)

% plot(independent_model.race_model_mean(independent_model.is_adjacent==1),...
%     independent_model.multi_mean_rt(independent_model.is_adjacent==1),'.','markersize',20,'color','k')
% hold on
% plot(independent_model.race_model_mean(independent_model.is_adjacent==0),...
%     independent_model.multi_mean_rt(independent_model.is_adjacent==0),'o','markersize',10,'color','k')

% adjacent
errorbar(independent_model.race_model_mean(independent_model.is_adjacent==1),...
    independent_model.multi_mean_rt(independent_model.is_adjacent==1),...
    independent_model.multi_std_err_rt(independent_model.is_adjacent==1),'vertical','.','markersize',20,'color','k')
errorbar(independent_model.race_model_mean(independent_model.is_adjacent==1),...
    independent_model.multi_mean_rt(independent_model.is_adjacent==1),...
    independent_model.race_model_std(independent_model.is_adjacent==1)./sqrt(num_trials),'horizontal','.','markersize',20,'color','k')

% far
errorbar(independent_model.race_model_mean(independent_model.is_adjacent==0),...
    independent_model.multi_mean_rt(independent_model.is_adjacent==0),...
    independent_model.multi_std_err_rt(independent_model.is_adjacent==0),'vertical','o','markersize',10,'color','k')
errorbar(independent_model.race_model_mean(independent_model.is_adjacent==0),...
    independent_model.multi_mean_rt(independent_model.is_adjacent==0),...
    independent_model.race_model_std(independent_model.is_adjacent==0)./sqrt(num_trials),'horizontal','o','markersize',10,'color','k')


hold on
plot([0,0.3],[0,0.3],'r--','linewidth',1.5)
xlabel('Race model RT (s)')
ylabel('Actual RT(s)')
formatForLee(gcf)
set(gca,'fontsize',16)
xlim([0.1,0.25])
ylim([0.1,0.25])
l = legend('adjacent','non-adjacent');
set(l,'box','off','location','best');


%% compare best electrode in a group to group RT

f=figure();
f.Name = [monkey_name,'_bestElecVsGroup'];

[min_single,min_idx] = cellfun(@min,independent_model.single_mean_rt);
min_single_std_err = [];
for i = 1:numel(independent_model.single_std_rt)
    min_single_std_err(end+1) = independent_model.single_std_rt{i}(min_idx(i))/sqrt(independent_model.single_num_trials{i}(min_idx(i)));
end

% legend placeholders
plot(-100,-100,'k.','markersize',20)
hold on
plot(-100,-100,'ko','markersize',10)

% adjacent
errorbar(min_single(independent_model.num_chans==2),...
    independent_model.multi_mean_rt(independent_model.num_chans==2),...
    min_single_std_err(independent_model.num_chans==2),'horizontal','.','markersize',20,'color','k')
    
errorbar(min_single(independent_model.num_chans==2),...
    independent_model.multi_mean_rt(independent_model.num_chans==2),...
    independent_model.multi_std_err_rt(independent_model.num_chans==2),'vertical','.','markersize',20,'color','k')

% far
errorbar(min_single(independent_model.num_chans==3),...
    independent_model.multi_mean_rt(independent_model.num_chans==3),...
    min_single_std_err(independent_model.num_chans==3),'horizontal','o','markersize',10,'color','k')
    
errorbar(min_single(independent_model.num_chans==3),...
    independent_model.multi_mean_rt(independent_model.num_chans==3),...
    independent_model.multi_std_err_rt(independent_model.num_chans==3),'vertical','o','markersize',10,'color','k')


plot([0,0.3],[0,0.3],'r--','linewidth',1.5)
xlabel('Fastest electrode RT (s)')
ylabel('Actual RT(s)')
formatForLee(gcf)
set(gca,'fontsize',16)
xlim([0.1,0.4])
ylim([0.1,0.4])
l = legend('Pair','Triplet');
set(l,'box','off','location','best');

%% compare best electrode in a group to race model RT
    f=figure();
f.Name = [monkey_name,'_bestElecVsRaceModel'];

plot(cellfun(@min,independent_model.single_mean_rt(independent_model.is_adjacent==1)),...
    independent_model.race_model_mean(independent_model.is_adjacent==1),'.','markersize',20,'color','k')
hold on
plot(cellfun(@min,independent_model.single_mean_rt(independent_model.is_adjacent==0)),...
    independent_model.race_model_mean(independent_model.is_adjacent==0),'o','markersize',10,'color','k')

hold on
plot([0,0.3],[0,0.3],'r--','linewidth',1.5)
xlabel('Fastest electrode RT (s)')
ylabel('Race Model RT(s)')
formatForLee(gcf)
set(gca,'fontsize',16)
xlim([0.1,0.25])
ylim([0.1,0.25])
l = legend('adjacent','non-adjacent');
set(l,'box','off','location','best');