function [ f ] = plotPsychometricCurve( psych_data_all,input_data )

    % check for enough colors
    if(numel(input_data.colors) < max(input_data.psych_data_idx_list))
        error('not enough colors')
    end
    
    % seperate bootstrap data from actual data
    psych_data = psych_data_all{1,input_data.axis};
    psych_data_boot = [];
    if(size(psych_data_all,1) > 1)
        psych_data_boot = psych_data_all(2:end,input_data.axis);
    end
    
    f=figure();
    hold on
    for i = input_data.psych_data_idx_list
        % make psychometric curve
        plot(psych_data(i).bump_dirs,psych_data(i).psych_curve_data,'.','color',input_data.colors{i},'markersize',14)
        xlabel('Bump Direction')
        ylabel('Percent to 0-deg target')
        xlim([0,180])
        ylim([-0.1,1.1])
        formatForLee(gcf)
        ax = gca;
        ax.FontSize = 14;
        ax.XTick = 0:45:180;
        ax.XAxis.MinorTickValues = [0:15:180];

        if(~isempty(psych_data(i).psych_fit))
            x_fit = [0:180];
            y_fit = feval(psych_data(i).psych_fit.fitObj,x_fit);
    %         y_fit_conf = predint(psych_data(i).psych_fit.fitObj,x_fit,0.95,'functional','off');

            hold on
            plot(x_fit,y_fit,'color',input_data.colors{i},'linewidth',1.5)
        end

        % if bootstrap data exists, add conf bounds
        if(~isempty(psych_data_boot) && input_data.plot_bootstrap)
            x_fit = [0:180];
            y_fit = zeros(numel(psych_data_boot),numel(x_fit));
            for boot = 1:numel(psych_data_boot)
                y_fit(boot,:) = feval(psych_data_boot{boot}(i).psych_fit.fitObj,x_fit)';
            end
            
            y_fit_sort = sort(y_fit,1);
            
            idx_use_upper = ceil(0.95*size(y_fit,1));
            idx_use_lower = floor(0.05*size(y_fit,1));
            
            y_fit_upper = y_fit_sort(idx_use_upper,:);
            y_fit_lower = y_fit_sort(idx_use_lower,:);
            
            hold on
            plot(x_fit,y_fit_lower,'color',input_data.colors{i},'linewidth',1.5,'linestyle','--')
            plot(x_fit,y_fit_upper,'color',input_data.colors{i},'linewidth',1.5,'linestyle','--')

        end
        
    end
  

end

