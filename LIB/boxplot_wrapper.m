function [  ] = boxplot_wrapper(x,data,input_data)
    % makes a single boxplot of data at position x
    
    % params are listed in setupParams
    
    % removes nan's
    data = data(~isnan(data));

    plot_data.median = median(data);
    plot_data.inter_quart = [prctile(data,25),prctile(data,75)];
    plot_data.outliers = data(data > plot_data.inter_quart(2) + 1.5*diff(plot_data.inter_quart) | data < plot_data.inter_quart(1) - 1.5*diff(plot_data.inter_quart));
    plot_data.whisker_pos = [min(data(data >= plot_data.inter_quart(1) - 1.5*diff(plot_data.inter_quart))),...
        max(data(data <= plot_data.inter_quart(2) + 1.5*diff(plot_data.inter_quart)))];

    plot_data.whisker_pos(1) = min(plot_data.whisker_pos(1),plot_data.inter_quart(1));
    plot_data.whisker_pos(2) = max(plot_data.whisker_pos(2),plot_data.inter_quart(2));
    
    % plot boxplot
    params = setupParams(input_data);
    hold on 
    % plot outliers
    if(~isempty(plot_data.outliers))
        plot(x,plot_data.outliers,'color',params.outlier_color,'marker',params.outlier_marker,'markersize',params.outlier_marker_size);
    end
    
    if(~params.use_log_x_scale)
        

        % plot box
        plot(x+params.box_width.*[-0.5,0.5],plot_data.inter_quart(2) + [0,0],'color', params.box_color,'linewidth',params.linewidth,'linestyle',params.box_linestyle);
        plot(x+params.box_width.*[-0.5,0.5],plot_data.inter_quart(1) + [0,0],'color', params.box_color,'linewidth',params.linewidth,'linestyle',params.box_linestyle);
        plot(x+params.box_width.*[-0.5,-0.5],[plot_data.inter_quart(1),plot_data.inter_quart(2)],'color', params.box_color,'linewidth',params.linewidth,'linestyle',params.box_linestyle);
        plot(x+params.box_width.*[0.5,0.5],[plot_data.inter_quart(1),plot_data.inter_quart(2)],'color', params.box_color,'linewidth',params.linewidth,'linestyle',params.box_linestyle);

        % plot whiskers
        plot(x+[0,0],[plot_data.inter_quart(2),plot_data.whisker_pos(2)],'color',params.whisker_color,'linewidth',params.linewidth,'linestyle',params.whisker_linestyle);
        plot(x+[0,0],[plot_data.inter_quart(1),plot_data.whisker_pos(1)],'color',params.whisker_color,'linewidth',params.linewidth,'linestyle',params.whisker_linestyle);
        plot(x+params.whisker_width*[-0.5,0.5],plot_data.whisker_pos(2)+[0,0],'color',params.whisker_color,'linewidth',params.linewidth,'linestyle',params.whisker_linestyle);
        plot(x+params.whisker_width*[-0.5,0.5],plot_data.whisker_pos(1)+[0,0],'color',params.whisker_color,'linewidth',params.linewidth,'linestyle',params.whisker_linestyle);

        % plot median line
        plot(x+params.box_width.*[-0.5,0.5],plot_data.median + [0,0],'color', params.median_color,'linewidth',params.linewidth,'linestyle',params.median_linestyle);
    elseif(params.use_log_x_scale)
       

        % plot box
        plot(x*[1./params.box_width,params.box_width],plot_data.inter_quart(2) + [0,0],'color', params.box_color,'linewidth',params.linewidth,'linestyle',params.box_linestyle);
        plot(x*[1./params.box_width,params.box_width],plot_data.inter_quart(1) + [0,0],'color', params.box_color,'linewidth',params.linewidth,'linestyle',params.box_linestyle);
        plot(x*[1./params.box_width,1./params.box_width],[plot_data.inter_quart(1),plot_data.inter_quart(2)],'color', params.box_color,'linewidth',params.linewidth,'linestyle',params.box_linestyle);
        plot(x*[params.box_width,params.box_width],[plot_data.inter_quart(1),plot_data.inter_quart(2)],'color', params.box_color,'linewidth',params.linewidth,'linestyle',params.box_linestyle);

        % plot whiskers
        plot(x+[0,0],[plot_data.inter_quart(2),plot_data.whisker_pos(2)],'color',params.whisker_color,'linewidth',params.linewidth,'linestyle',params.whisker_linestyle);
        plot(x+[0,0],[plot_data.inter_quart(1),plot_data.whisker_pos(1)],'color',params.whisker_color,'linewidth',params.linewidth,'linestyle',params.whisker_linestyle);
        plot(x*[1./params.whisker_width,params.whisker_width],plot_data.whisker_pos(2)+[0,0],'color',params.whisker_color,'linewidth',params.linewidth,'linestyle',params.whisker_linestyle);
        plot(x*[1./params.whisker_width,params.whisker_width],plot_data.whisker_pos(1)+[0,0],'color',params.whisker_color,'linewidth',params.linewidth,'linestyle',params.whisker_linestyle);     
        
        % plot median line
        plot(x*[1./params.box_width,params.box_width],plot_data.median + [0,0],'color', params.median_color,'linewidth',params.linewidth,'linestyle',params.median_linestyle);
    end
    
end

function [params] = setupParams(input_params)

    params.linewidth = 2;
    params.box_width = 0.5;
    
    params.use_same_color_for_all = 0;
    params.master_color = 'k';
    
    params.median_linestyle = '-';
    params.median_color = 'k';
    
    params.box_color = 'k';
    params.box_linestyle = '-';
    
    params.whisker_color = 'k';
    params.whisker_linestyle = '-';
    params.whisker_width = 0.2;
    
    params.outlier_marker = 'x';
    params.outlier_marker_size = 8;
    params.outlier_color = 'k';
    
    params.use_log_x_scale = 0;
    %% check if in optsSave and optsSaveInput, overwrite if so
    try
        inputFieldnames = fieldnames(input_params);
        for fn = 1:numel(inputFieldnames)
           if(isfield(input_params,inputFieldnames{fn}))
               params.(inputFieldnames{fn}) = input_params.(inputFieldnames{fn});
           end
        end
    catch
        % do nothing, [] was inputted which means use default setting
    end
    
    if(params.use_same_color_for_all == 1)
        
        params.box_color = params.master_color;
        params.whisker_color = params.master_color;
        params.outlier_color = params.master_color;
    end
    
    
end