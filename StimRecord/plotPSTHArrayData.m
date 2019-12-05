function [output_data] = plotPSTHArrayData(array_data,input_data)

    %plots psth's for each condition in array_data
    max_y_lim = 0;
    f_list = {};

    for u = input_data.unit_idx
        max_y_lim = 0;
        f_list{u} = figure();
        f_list{u}.Position = [251.4000 -49.4000 704 500*input_data.num_cols];
        f_list{u}.Name = [array_data{u}.monkey,'_DblPulseTrains_PSTH_chan',num2str(array_data{u}.CHAN_REC),'_',num2str(array_data{u}.ID),'ID'];%'rec_IPI',num2str(input_data.IPI(cond)),...
%                 '_numPulses',num2str(input_data.num_pulses(cond)),'_unitID',num2str(u)];

        if(isfield(input_data,'suffix'))
            f_list{u}.Name = [f_list{u}.Name,'_',input_data.suffix];
        end
        for condition = 1:numel(array_data{u}.binCounts)
            subplot(ceil(numel(array_data{u}.binCounts)/input_data.num_cols),input_data.num_cols,condition)
%             subplot(ceil(numel(array_data{u}.binCounts)/input_data.num_cols),1,ceil(condition/input_data.num_cols))
            
            if(exist('downsample_stims') && downsample_stims == 1) % use num_stims_use 
%                 plot(array_data{u}.binEdges{cond}(1:end-1)+mean(diff(array_data{u}.binEdges{cond}))/2,...
%                     array_data{u}.binCounts{cond}/(mode(diff(array_data{u}.binEdges{cond}/1000)))/min(num_stims_use,array_data{u}.num_stims(cond)),'color',getColorFromList(1,1),'linewidth',1.5)
            else % use num_stims from array_data
                plot(array_data{u}.binEdges{condition}(1:end-1)+mean(diff(array_data{u}.binEdges{condition}))/2,...
                    array_data{u}.binCounts{condition}/(mode(diff(array_data{u}.binEdges{condition}/1000))),...
                        'color',getColorFromList(1,1),'linewidth',1.5);
%                         'color',getColorFromList(1,mod(condition,input_data.num_cols)),'linewidth',1.5)
                
                max_y_lim = max(max_y_lim,max(array_data{u}.binCounts{condition}/(mode(diff(array_data{u}.binEdges{condition}/1000)))));
                hold on
                stim_on_line = plot([0,0],[0,1000],'r--');
                stim_off_line = plot([4000,4000],[0,1000],'r--');
                uistack(stim_on_line,'bottom');
                uistack(stim_off_line,'bottom');
            end
            formatForLee(gcf)
            if(condition == 1)
                xlabel('Time after 1st pulse (ms)');
                ylabel('Firing rate (spks/s)')
            end
            set(gca,'fontsize',14);
            
        end
        
        for condition = 1:numel(f_list{u}.Children)
            f_list{u}.Children(condition).YLim(end) = max_y_lim;
            f_list{u}.Children(condition).XLim = input_data.window;
        end
    end



end