function [plot_handle] = errorbar_LIB(x_vals,y_means,y_std_err,varargin)

    % set default parameters
    marker = '.';
    markersize = 16;
    linewidth = 1.5;
    color = 'k';
    linestyle = '-';
    caplength = 0.3; % 0.3 of an x unit
    
    % override default if in varargin
    for i = 1:2:numel(varargin)
        switch varargin{i}
            case 'marker'
                marker = varargin{i+1};
            case 'markersize'
                markersize = varargin{i+1};
            case 'linewidth'
                linewidth = varargin{i+1};
            case 'color'
                color = varargin{i+1};
            case 'linestyle'
                linestyle = varargin{i+1};
            case 'caplength'
                caplength = varargin{i+1};
        end
    end

    % plot line for error bar
    hold on
    for i = 1:numel(x_vals)
        % vertical line
        plot([x_vals(i),x_vals(i)],y_means(i)+y_std_err(i)*[-1,1],'linewidth',linewidth,'linestyle',linestyle,'color',color);
        
        % cap
        plot(x_vals(i)+caplength*[-1,1],y_means(i)+y_std_err(i)*[-1,-1],'linewidth',linewidth,'linestyle',linestyle,'color',color);
        plot(x_vals(i)+caplength*[-1,1],y_means(i)+y_std_err(i)*[1,1],'linewidth',linewidth,'linestyle',linestyle,'color',color);
        
    end
    % plot marker for means
    plot(x_vals,y_means,'marker',marker,'markersize',markersize,'color',color);
    
    
    % done




end