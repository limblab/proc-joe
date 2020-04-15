function [output_data] = plotPSTHArrayData(array_data,input_data)

    %plots psth's for each condition in array_data
    max_y_lim = 0;
    f_list = {};
%     amp_freq_plot_map = [10,11,12,7,8,9,4,5,6,1,2,3];
%     intermittent_plot_map = [3,2,1,6,5,4,9,8,7,12,11,10];
    amp_freq_plot_map = [1,2,3,4,5,6,7,8];
    for u = input_data.unit_idx
        max_y_lim = 0;
        f_list{u} = figure();

        f_list{u}.Position = [-863 -61.4000 707.2000 788];
%         f_list{u}.Name = [array_data{u}.monkey,'_longTrains_PSTH_chan',num2str(array_data{u}.CHAN_SENT{1}),'stim_chan'...
%             num2str(array_data{u}.CHAN_REC),'rec_',num2str(array_data{u}.ID),'ID'];%'rec_IPI',num2str(input_data.IPI(cond)),...
% %                 '_numPulses',num2str(input_data.num_pulses(cond)),'_unitID',num2str(u)];

        if(isfield(input_data,'suffix'))
            f_list{u}.Name = [f_list{u}.Name,'_',input_data.suffix];
        end
        for condition = 1:numel(array_data{u}.binCounts)
            subplot_idx = condition;
%             if(input_data.amp_freq)
%                 subplot_idx = amp_freq_plot_map(condition);
%             else
%                 subplot_idx = intermittent_plot_map(condition);
%             end
            subplot(ceil(numel(array_data{u}.binCounts)/input_data.num_cols),input_data.num_cols,subplot_idx)
%             subplot(ceil(numel(array_data{u}.binCounts)/input_data.num_cols),1,ceil(condition/input_data.num_cols))
            
            if(exist('downsample_stims') && downsample_stims == 1) % use num_stims_use 
%                 plot(array_data{u}.binEdges{cond}(1:end-1)+mean(diff(array_data{u}.binEdges{cond}))/2,...
%                     array_data{u}.binCounts{cond}/(mode(diff(array_data{u}.binEdges{cond}/1000)))/min(num_stims_use,array_data{u}.num_stims(cond)),'color',getColorFromList(1,1),'linewidth',1.5)
            else % use num_stims from array_data
                x_data = array_data{u}.binEdges{condition}(1:end-1)+mean(diff(array_data{u}.binEdges{condition}))/2;
                % if bin counts are in number of spikes (instead of mean
                % number of spikes), divide by num stims
                if(sum(floor(array_data{u}.binCounts{condition}) ~= array_data{u}.binCounts{condition}) == 0)
                    array_data{u}.binCounts{condition} = array_data{u}.binCounts{condition}/array_data{u}.num_stims(condition);
                end
                if(input_data.account_for_artifact)
                    num_stims_in_bin = histcounts(array_data{u}.PULSE_TIMES{1}{1}*1000,array_data{u}.binEdges{condition});
                    
                    y_data = array_data{u}.binCounts{condition}./(mode(diff(array_data{u}.binEdges{condition}/1000))-num_stims_in_bin/1000); % remove 1 ms per stim
                else
                    y_data = array_data{u}.binCounts{condition}/(mode(diff(array_data{u}.binEdges{condition}/1000)));
                end
                plot(x_data,y_data,...
                        'color',getColorFromList(1,1),'linewidth',1.5);
                
                max_y_lim = max(max_y_lim,max(y_data));
                hold on
                if(isfield(array_data{u},'PULSE_TIMES') && ~isempty(array_data{u}.PULSE_TIMES{condition}))
                    stim_on_line = plot([0,0],[0,1000],'r--');
                    stim_off_line = plot([max(array_data{u}.PULSE_TIMES{condition}{1})*1000,max(array_data{u}.PULSE_TIMES{condition}{1})*1000],[0,1000],'r--');
                    uistack(stim_on_line,'bottom');
                    uistack(stim_off_line,'bottom');
                end
            end
            formatForLee(gcf)
            if(subplot_idx == 1)
                xlabel('Time after 1st pulse (s)');
                ylabel('Firing rate (spks/s)')
            end
            set(gca,'fontsize',14);
            ax=gca;
%             ax.XTick = [-4000,0,4000,8000,12000];
%             ax.XTickLabel = {'-4','0','4','8','12'};
            ax.XMinorTick = 'on';
        end
        
        for condition = 1:numel(f_list{u}.Children)
            f_list{u}.Children(condition).YLim(end) = max_y_lim;
            f_list{u}.Children(condition).XLim = input_data.window;
        end
    end



end