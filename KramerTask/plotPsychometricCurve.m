function [ f ] = plotPsychometricCurve( psych_data,input_data )

    if(numel(input_data.colors) < max(input_data.psych_data_idx_list))
        error('not enough colors')
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
%         plot(x_fit,y_fit_conf,'color',input_data.colors{i},'linewidth',1.5,'linestyle','--');
        
%         if(isfield(psych_data(i),'bootstrap_fit_params')) % deal with the bootstrapped data
%             % find 2.5 percentile, 97.5 percentile parameters
%             [sort_param_list] = sort(psych_data(i).bootstrap_fit_params);
%             conf_idx = ceil(size(psych_data(i).bootstrap_fit_params,1)*[0.025,0.975]);
%             conf_params = psych_data(i).bootstrap_fit_params(conf_idx,:);
%             
%             x_fit = [0:180];
%             y_fit = conf_params(:,1) + conf_params(:,2).*(erf(conf_params(:,3).*(x_fit-conf_params(:,4))));
%             
%             plot(x_fit,y_fit,'color',input_data.colors{i},'linewidth',1.5,'linestyle','--')
%             
%         end
    end
  

end

