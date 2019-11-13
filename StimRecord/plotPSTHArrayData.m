function [output_data] = plotPSTHArrayData(array_data,input_data)

    %plots psth's for each condition in array_data
    max_y_lim = 0;
    f_list = {};

    for u = input_data.unit_idx
        max_y_lim = 0;
        f_list{u} = figure();
        f_list{u}.Position = [251.4000 -49.4000 704 804.8000];
        f_list{u}.Name = [array_data{u}.monkey,'_DblPulseTrains_PSTH_chan',num2str(array_data{u}.CHAN_REC)];%'rec_IPI',num2str(input_data.IPI(cond)),...
%                 '_numPulses',num2str(input_data.num_pulses(cond)),'_unitID',num2str(u)];

        if(isfield(input_data,'suffix'))
            f_list{u}.Name = [f_list{u}.Name,'_',input_data.suffix];
        end
        for cond = 1:numel(array_data{u}.binCounts)
            subplot(ceil(numel(array_data{u}.binCounts)/2),2,cond)

            
            if(exist('downsample_stims') && downsample_stims == 1) % use num_stims_use 
                plot(array_data{u}.binEdges{cond}(1:end-1)+mean(diff(array_data{u}.binEdges{cond}))/2,...
                    array_data{u}.binCounts{cond}/(mode(diff(array_data{u}.binEdges{cond}/1000)))/min(num_stims_use,array_data{u}.num_stims(cond)),'color',getColorFromList(1,1),'linewidth',1.5)
            else % use num_stims from array_data
                plot(array_data{u}.binEdges{cond}(1:end-1)+mean(diff(array_data{u}.binEdges{cond}))/2,...
                    array_data{u}.binCounts{cond}/(mode(diff(array_data{u}.binEdges{cond}/1000)))/array_data{u}.num_stims(cond),'color',getColorFromList(1,1),'linewidth',1.5)
            end
            formatForLee(gcf)
            xlabel('Time after 1st pulse (ms)');
            ylabel('Firing rate (spks/s)')
            set(gca,'fontsize',14);
            max_y_lim = max(max_y_lim,f_list{u}.Children(1).YLim(end));
        end
        
        for cond = 1:numel(array_data{u}.binCounts)
            f_list{u}.Children(cond).YLim(end) = max_y_lim;
            f_list{u}.Children(cond).XLim = input_data.window;
        end
    end



end