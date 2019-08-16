function [output_data] = plotPSTHArrayData(array_data,input_data)

    %plots psth's for each condition in array_data
    max_y_lim = 0;
    for u = 1:numel(input_data.unit_idx)
        max_y_lim = 0;
        f_list = {};
        for cond = 1:numel(array_data{u}.binCounts)
            f_list{cond} = figure();
            f_list{cond}.Name = ['Han_20190401_chan',num2str(input_data.chan_rec),'rec_IPI',num2str(input_data.IPI(cond)),...
                '_numPulses',num2str(input_data.num_pulses(cond)),'_unitID',num2str(u)];
            if(isfield(input_data,'suffix'))
                f_list{cond}.Name = [f_list{cond}.Name,'_',input_data.suffix];
            end
            
            if(exist('downsample_stims') && downsample_stims == 1) % use num_stims_use 
                plot(array_data{u}.binEdges{cond}(1:end-1)+mean(diff(array_data{u}.binEdges{cond}))/2,...
                    array_data{u}.binCounts{cond}/(input_data.bin_size/1000)/min(num_stims_use,array_data{u}.num_stims(cond)))
            else % use num_stims from array_data
                plot(array_data{u}.binEdges{cond}(1:end-1)+mean(diff(array_data{u}.binEdges{cond}))/2,...
                    array_data{u}.binCounts{cond}/(input_data.bin_size/1000)/array_data{u}.num_stims(cond))
            end
            formatForLee(gcf)
            xlabel('Time after 1st pulse (ms)');
            ylabel('Firing rate (spks/s)')
            set(gca,'fontsize',14);
            max_y_lim = max(max_y_lim,f_list{cond}.Children.YLim(end));
        end
        
        for cond = 1:numel(array_data{u}.binCounts)
            f_list{cond}.Children.YLim(end) = max_y_lim;
            f_list{cond}.Children.XLim = input_data.window;
        end
    end



end