%% get median and IQR for each condition, make plot
    amp_offset = [-0.15,-0.15,0.15,0.15]*4;
    
    amp_mask = mod(amps,10) == 0;
    
    median_list = zeros(sum(amp_mask),4);
    prctile_list = zeros(sum(amp_mask),4,2);
    
    f = figure();
    counter = 0;
    for a = 1:numel(amps)
        if(amp_mask(a) == 1)
            counter = counter + 1;
            % BR cath
            median_list(counter,1) = median(time_off_rails_cathodic(a,~isnan(time_off_rails_cathodic(a,:))));
            prctile_list(counter,1,:) = prctile(time_off_rails_cathodic(a,~isnan(time_off_rails_cathodic(a,:))),[25,75]);

            % duke cath
            median_list(counter,2) = median(time_recovered{a,1});
            prctile_list(counter,2,:) = prctile(time_recovered{a,1},[25,75]);

            % BR anod
            median_list(counter,3) = median(time_off_rails_anodic(a,~isnan(time_off_rails_anodic(a,:))));
            prctile_list(counter,3,:) = prctile(time_off_rails_anodic(a,~isnan(time_off_rails_anodic(a,:))),[25,75]);

            % duke anod
            median_list(counter,4) = median(time_recovered{a,2});
            prctile_list(counter,4,:) = prctile(time_recovered{a,2},[25,75]);
        end
    end 

    
    markers = {'s','.','s','.'};
    marker_sizes = [8,22,8,20];
    line_styles = {'--','-','--','-'};
    
    colors = [getColorFromList(1,0); getColorFromList(1,0); getColorFromList(1,1); getColorFromList(1,1)];
    
    for cond = 1:4
        errorbar(amps(amp_mask)+amp_offset(cond),median_list(:,cond),...
            median_list(:,cond) - squeeze(prctile_list(:,cond,1)),squeeze(prctile_list(:,cond,2)-median_list(:,cond)),...
            'linewidth',1.5,'marker',markers{cond},'markersize',marker_sizes(cond),'color',colors(cond,:),'linestyle',line_styles{cond})
        hold on
    end


    ylim([0,4]);
    xlim([0,105]);
    formatForLee(gcf);
    xlabel('Amplitude (\muA)');
    ylabel('Recovery time (ms)');
    ax = gca;
    set(ax,'fontsize',14)
    ax.XTick = [0:20:100];
    ax.XAxis.MinorTickValues = [0:10:100];













