%% setup code maybe

    td_all = td_all(103:end);
    td_bump_block = td_all(1:69);
    td_vis_block = td_all(70:end);
    td_vis = td_vis_block([td_vis_block.isBumpTrial] == 1);
    td_bump = td_bump_block;

%% given td_vis and td_bump, get mean firing rates during bump for each unit
    window = [0,0.2]; % s
    baseline = [-0.3,-0.1]; % s
    num_units = size(td_bump(1).LeftS1_unit_guide,1);
    total_spikes = [];
    total_spikes.vis = zeros(num_units,1,3);
    total_spikes.bump = zeros(num_units,1,3);
    total_spikes.num_trials = zeros(2,1); % baseline computed from bump block?
    
    for i = 1:2 % vis block = 1, bump block = 2
        switch i
            case 1
                td_use = td_vis;
                field_name = 'vis';
            case 2 
                td_use = td_bump;
                field_name = 'bump';
        end
        for tr = 1:numel(td_use)
            if(~isnan(td_use(tr).idx_goCueTime))
                window_idx = td_use(tr).idx_goCueTime + floor(window/td_use(tr).bin_size);
                baseline_idx = td_use(tr).idx_goCueTime + floor(baseline/td_use(tr).bin_size);
                total_spikes.(field_name)(:,total_spikes.num_trials(i)+1,i) = sum(td_use(tr).LeftS1_spikes(window_idx(1):window_idx(2),:))';
                total_spikes.(field_name)(:,total_spikes.num_trials(i)+1,3) = sum(td_use(tr).LeftS1_spikes(baseline_idx(1):baseline_idx(2),:))';
                total_spikes.num_trials(i) = total_spikes.num_trials(i) + 1;

            end
        end
            
        total_spikes.mean_fr(:,i) = mean(squeeze(total_spikes.(field_name)(:,:,i))')/(diff(window));
        total_spikes.std_fr(:,i) = std(squeeze(total_spikes.(field_name)(:,:,i))'/(diff(window)));
        if(i==2)
            total_spikes.mean_fr_baseline(:,1) = mean(squeeze(total_spikes.(field_name)(:,:,3))')/(diff(baseline));
            total_spikes.std_fr_baseline(:,1) = std(squeeze(total_spikes.(field_name)(:,:,3))'/diff(baseline));
        end
    end
    
%% plot each neuron's mean fr during the bump in each block
    plot_mask = total_spikes.mean_fr(:,2) > total_spikes.mean_fr_baseline + 0.5*total_spikes.std_fr_baseline;
    figure();
    plot(mean_fr(plot_mask,1),mean_fr(plot_mask,2),'k.','markersize',12)
    hold on
    plot([0,40],[0,40],'r--','linewidth',2)
    xlabel('Vis block (spk/s)')
    ylabel('Bump block (spk/s)')
    formatForLee(gcf)
    set(gca,'fontsize',14)
    
   