function [figs] = plotBinnedSpikeTimes(mdl_data, exp_data, input_data)

    % plots data comparing spike times for mdl vs exp
    figs = [];

    % latency of activation -- plot sum of spikes in each bin -- experiment
    % vs each of 3 diamters, across cell types
    
    figs(end+1) = figure('Position',[2097 219 1552 714]); hold on;
    color_list = inferno(numel(input_data.amps_plot) + 2); % + 2 to remove yellows
    ax_list = [];
    
    for data_type = 1:4 % experiment and 3 diameters
        if(data_type == 1) % experiment
            data = exp_data;
        else % model
            data = mdl_data;
        end
        for i_amp = 1:numel(input_data.amps_plot)
            amp_idx = find([data.amp] == input_data.amps_plot(i_amp),1,'first');
            
            % get mask (diameter)
            if(data_type == 1) % experiment
                mask = ones(size(data(amp_idx).binned_data,1),1);
            else
                mask = data(amp_idx).diam_list == data_type-1;
            end
            
            if(~isempty(amp_idx))
                ax_list(end+1,1) = subplot(2,4,data_type); hold on;
                plot(data(amp_idx).bin_centers*1000,mean(data(amp_idx).binned_data(mask==1,:)),'color',color_list(i_amp,:),'linewidth',2);
            end
        end
        
        % horiz line at 0
        plot([exp_data(1).bin_edges(1),exp_data(1).bin_edges(end)]*1000,[0,0],'k--')
        % shade blank period (before 1 ms)
        patch([-0.05,1,1,-0.05,-0.05],[-0.05,-0.05,1,1,-0.05],'r','FaceColor',getColorFromList(1,0),'FaceAlpha',0.25,'EdgeColor','none');

        formatForLee(gcf);
        set(gca,'fontsize',14)
        xlabel('Time after stim offset (ms)');
        if(data_type ==1)
            ylabel('Num spikes per stim per neuron');
        end
    end
        
    linkaxes(ax_list,'xy');
    xlim([exp_data(1).bin_edges(1),exp_data(1).bin_edges(end)]*1000);
    ylim([-0.05,1]);


    % plot number of spikes in first 2 ms bins for the model -- across 3
    % diameters, plot each cell type
    ax_list = [];
    % experiment. Plot data within specified windows -- first window is a
    % solid line
    linestyle_list = {'-','--','-.',':'};
    
    ax_list(1) = subplot(2,4,5); hold on;
    amps_plot = [exp_data.amp];
    exp_prob_spike = zeros(size(input_data.exp_windows,1),numel(amps_plot));
    for i_wind = 1:size(input_data.exp_windows,1)
        for i_amp = 1:numel(amps_plot)
            amp_idx = find([exp_data.amp] == amps_plot(i_amp),1,'first');
            bin_idx = find(exp_data(amp_idx).bin_centers > input_data.exp_windows(i_wind,1) & exp_data(amp_idx).bin_centers <= input_data.exp_windows(i_wind,2));
            exp_prob_spike(i_wind,i_amp) = mean(sum(exp_data(amp_idx).binned_data(:,bin_idx),2));
        end
        
        plot(amps_plot,exp_prob_spike(i_wind,:),'color','k','linestyle',linestyle_list{i_wind},'linewidth',2)
    end
    formatForLee(gcf);
    set(gca,'fontsize',14)
    xlabel('Amplitude (\muA)');
    ylabel('Num spikes per stim per neuron');
    
    
    % model
    amps_plot = [mdl_data.amp];
    mdl_prob_spike = zeros(3,numel(input_data.cell_ids),numel(amps_plot)); % 3 diameters
    for data_type = 1:3 % 3 diameters
        ax_list(data_type+1) = subplot(2,4,data_type+5); hold on;
        for i_cell_id = 1:numel(input_data.cell_ids)
            for i_amp = 1:numel(amps_plot)
                amp_idx = find([mdl_data.amp] == amps_plot(i_amp),1,'first');
                % get mask (diameter) and cell id
                mask = mdl_data(amp_idx).diam_list == data_type & mdl_data(amp_idx).cell_id_list == input_data.cell_ids(i_cell_id);
                bin_idx = find(mdl_data(amp_idx).bin_centers <= input_data.max_mdl_time);
                if(~isempty(amp_idx))
                    mdl_prob_spike(data_type,i_cell_id,i_amp) = mean(sum(mdl_data(amp_idx).binned_data(mask==1,bin_idx),2));
                end
            end
        
            % plot data
            plot(amps_plot,squeeze(mdl_prob_spike(data_type,i_cell_id,:)),'color',getColorFromList(1,i_cell_id+2),'linewidth',2);
        end
        formatForLee(gcf);
        set(gca,'fontsize',14)
        xlabel('Amplitude (\muA)');
    end
        
    linkaxes(ax_list,'xy');
    xlim([0,100]);
    ylim([-0.05,1]);

end