%% set up code:chan mapping
filename = 'Han_20181126';

mapping = [0,26;...
            1,61;...
            2,44;...
            3,11;...
            4,3];
%% make struct with mean_rt, std_rt and chan num
    curr_file_data = [];
    curr_file_data.chan = [];
    curr_file_data.mean_rt = [];
    curr_file_data.std_rt = [];
    curr_file_data.filename = {};
    curr_file_data.num_trials = [];
    curr_file_data.rt_all = {};
    
    for i = 1:size(mapping,1)
        % find stim code and bump mag == 0
        cue_idx = find([data.cueInfo.bumpMag] == 0 & [data.cueInfo.stimCode] == mapping(i,1));
        if(~isempty(cue_idx))
            curr_file_data.chan(end+1) = mapping(i,2);
            curr_file_data.mean_rt(end+1) = mean(data.cueInfo(cue_idx).rt);
            curr_file_data.std_rt(end+1) = std(data.cueInfo(cue_idx).rt);
            curr_file_data.num_trials(end+1) = numel(data.cueInfo(cue_idx).rt);
            curr_file_data.filename{end+1} = filename;
            curr_file_data.rt_all{end+1} =  data.cueInfo(cue_idx).rt;
        end
    end
  
%% load in struct with all of those, and combine
    all_files_data.chan = [all_files_data.chan, curr_file_data.chan];
    all_files_data.mean_rt = [all_files_data.mean_rt, curr_file_data.mean_rt];
    all_files_data.std_rt = [all_files_data.std_rt, curr_file_data.std_rt];
    all_files_data.filename = [all_files_data.filename,curr_file_data.filename];
    all_files_data.num_trials = [all_files_data.num_trials,curr_file_data.num_trials];
    all_files_data.rt_all = [all_files_data.rt_all,curr_file_data.rt_all];
    all_files_data
    
    
%% dot plot of RTs (instead of a histogram?)
    f=figure();
    hold on;
    f.Position = [278.6000 95.4000 560 420];
    f.Name = 'Han_singleElectrodes_dotPlot';
    ax=gca;
    % get indexes of each day
    filename_list = unique(all_files_data.filename);
    end_idx = 0; % will be set to start idx in for loop, start at beginning
    for fnum = 1:numel(filename_list)
        % get data for this day
        start_idx = end_idx+1;
        end_idx = max(find(strcmpi(all_files_data.filename,filename_list{fnum})));
        x_data = (start_idx:1:end_idx) + fnum - 1;
        
        % plot bump bar for this day
        % bump data [0.1503,0.1401,0.1489,0.1450,0.1476,0.1659,0.1550,0.1639,0.1650,0.1603,0.1439] for [10/15,10/16,10/17,10/19,10/22,11/15,11/15,11/19,11/21,11/27,11/28] in Han
        % std err bump data [0.0017,0.003,0.0029,0.0032,0.0035,0.0028,0.0027,0.0022,0.0015,0.003,0.0034] for
        % [10/15,10/16,10/17,10/19,10/22,11/15,11/16,11/19,11/21,11/27] in Han
        bump_list = [0.1503,0.1401,0.1489,0.1450,0.1476,0.1659,0.1550,0.1639,0.1650,0.1603,0.1439]; % Han
        bump_std_list = [0.0017,0.003,0.0029,0.0032,0.0035,0.0028,0.0027,0.0022,0.0015,0.003,0.0034]; % Han
        bump_bar = fill([x_data(1)-1,x_data(end)+1,x_data(end)+1,x_data(1)-1,x_data(1)-1],...
            bump_list(fnum)+bump_std_list(fnum).*[-1,-1,1,1,-1],...
            'k','EdgeColor','none','FaceAlpha',0.5);
        
        uistack(bump_bar,'bottom')
        plot([x_data(1)-1,x_data(end)+1],[bump_list(fnum),bump_list(fnum)],'k','linewidth',1.5)

        % plot visual bar for this day
        % visual data [0.2141 0.2004 0.2170 0.2281, 0.2178,0.2199,0.235,0.2122,0.247,0.2212,0.2133] for [10/15,10/16,10/17,10/19,10/22,11/15,11/16,11/19,11/21,11/27,11/28] in Han
        % std err vis data [0.0076,0.0099,0.0074,0.0089,0.0109,0.0041,0.0105,0.0043,0.0068,0.0075,0.0074] for
        % [10/15,10/16,10/17,10/19,10/22,11/15,11/16,11/19,11/21,11/27,11/28] in Han
        vis_list = [0.2141 0.2004 0.2170 0.2281, 0.2178,0.2199,0.235 0.2122,0.247,0.2212,0.2133]; % Han
        vis_std_list = [0.0076,0.0099,0.0074,0.0089,0.0109,0.0041,0.0105 0.0043,0.0068,0.0075,0.0074]; % Han
        vis_bar = fill([x_data(1)-1,x_data(end)+1,x_data(end)+1,x_data(1)-1,x_data(1)-1],...
            vis_list(fnum)+vis_std_list(fnum).*[-1,-1,1,1,-1],...
            'k','EdgeColor','none','FaceAlpha',0.5);
        
        uistack(vis_bar,'bottom')
        plot([x_data(1)-1,x_data(end)+1],[vis_list(fnum),vis_list(fnum)],'k--','linewidth',1.5)

        
        % plot data for this day
        plot(x_data,all_files_data.mean_rt(start_idx:end_idx),'.','color',getColorFromList(1,0),'markersize',18) 
        for i = 1:numel(x_data)
            errorbar(x_data(i),all_files_data.mean_rt(i+start_idx-1),all_files_data.std_rt(i+start_idx-1)./sqrt(all_files_data.num_trials(i+start_idx-1)),'color',getColorFromList(1,0),'linewidth',1.5);
        end
        % plot dashed line to denote end of a day
        if(fnum < numel(filename_list)) % don't do this for the last one
            plot(x_data(end)+1+[0,0],[0,1],'--','color',[0.5,0.5,0.5],'linewidth',1.5)
        end
        
    end
    
    ylim([0,0.4])
    xlim([0,x_data(end)+1])
%     set(gca,'XTickLabel',{});
%     set(gca,'XTick',[]);
    xlabel('Electrode');
    ylabel('RT (s)')
    formatForLee(gcf)
    ax.FontSize = 14;
    
    
    
%% vertical histogram of RTs (meant to go next to above plot)
    binSize = 0.01;
    binEdges = 0.1:binSize:0.4;
    f=figure();
    f.Name = 'Han_singleElectrodes_vertHist';
    f.Position = [730.6000 155.4000 221.6000 420];
    ax=gca;
    binCounts = histcounts(all_files_data.mean_rt,binEdges);
    h = barh(binEdges(1:end-1)+binSize/2,binCounts,'BarWidth',1);
    h.FaceAlpha = 1;
    h.FaceColor = getColorFromList(1,0);
    h.EdgeColor = 'k';
    ax.YAxis.TickLabels = {};
    
    xlabel('Number of electrodes');
%     ylabel('RT (s)');
    formatForLee(gcf)
    set(gca,'fontsize',14)
    ylim([0.,0.4])

%% plot heatmap of RT across array
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
    map_data = loadMapFile(mapFileName);
    c_map = inferno();
    c_map = flip(c_map,1);
    rt_min = 0.1;
    rt_max = 0.35;
    f=figure();
    f.Name = 'Han_singleElectrodes_arrayMap';
    hold on
    for chan_idx = 1:numel(all_files_data.chan)
        map_idx = find(all_files_data.chan(chan_idx) == map_data.chan);
        mean_rt = all_files_data.mean_rt(chan_idx);
        c_map_idx = ceil((mean_rt-rt_min)/(rt_max-rt_min)*size(c_map,1));
        rt_color = c_map(min(max(c_map_idx,1),size(c_map,1)),:);
        rectangle('Position',[map_data.row(map_idx),11-map_data.col(map_idx),1,1],'FaceColor',rt_color, 'EdgeColor','none');
        hold on
    end
    plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5)

    % plot bad channels
    for i = 1:numel(all_files_data.no_response)
        map_idx = find(all_files_data.no_response(i) == map_data.chan);
        pos = [map_data.row(map_idx),11-map_data.col(map_idx)];
        plot([pos(1),pos(1)+1],[pos(2),pos(2)+1],'k-','linewidth',2)
        plot([pos(1),pos(1)+1],[pos(2)+1,pos(2)],'k-','linewidth',2)
    end
    
    f=gcf;
    set(gca,'visible','off')
    xlim([1,11])
    ylim([1,11])
    axis square
    b=colorbar();
    colormap(c_map);
    b.FontSize = 14;
    b.Ticks = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]; 
    b.TickDirection = 'out';
    b.TickLabels = cell(1,numel(b.Ticks));
    for i = 1:2:numel(b.Ticks)
        if(i==numel(b.Ticks))
            b.TickLabels{i,1} = strcat(num2str((i-1)*(rt_max-rt_min)/(numel(b.Ticks)-1) + rt_min));
        elseif(i==1)
            b.TickLabels{i,1} = strcat(num2str((i-1)*(rt_max-rt_min)/(numel(b.Ticks)-1) + rt_min));
        else
            b.TickLabels{i,1} = num2str((i-1)*(rt_max-rt_min)/(numel(b.Ticks)-1) + rt_min);
        end

    end
    b.Label.String = 'RT (s)';
    b.Label.FontSize = 16;
    
%% plot RT vs impedance of electrodes. 
%skip the first data as that is the 'internal channel' of the
    %stimulator. data for electrodes 1:96 are in rows 2:97 of the stim
    %impedance object.
    imp.impedance = imp.impedance(2:97)/1000;
    f=figure();
    f.Name = 'Han_singleElectrodes_impedance';
    plot(imp.impedance(all_files_data.chan),all_files_data.mean_rt,'k.','markersize',16)
    formatForLee(gcf)
    xlabel('Impedance (kOhm)');
    ylabel('RT (s)');
    set(gca,'fontsize',14)
            
            
%% plot mean rt vs std rt for each channel
    f = figure();
    f.Name = 'Han_singleElectrodes_meanVsStd';
    plot(all_files_data.mean_rt,all_files_data.std_rt,'k.','markersize',16)
    hold on
    [mean_std_fit,mean_std_gof] = fit(all_files_data.mean_rt',all_files_data.std_rt','a*x+b');
    plot([0.1,0.35],...
        feval(mean_std_fit,[0.1,0.35]),...
        'r--','linewidth',1.5)
    
    
    formatForLee(gcf);
    xlabel('Mean RT (s)');
    ylabel('Std RT (s)');
    set(gca,'fontsize',14)
    xlim([0,0.4]);
    ylim([0,0.15]);
    
%% plot histogram of std rt's
    f = figure();
    f.Name = 'Han_singleElectrodes_stdHist';
    bin_edges = 0:0.01:0.15;
    histogram(all_files_data.std_rt,bin_edges);
    formatForLee(gcf);
    xlabel('Std RT (s)');
    ylabel('Number of electrodes');
    set(gca,'fontsize',14);
            
%% plot slope histogram
    slopes = [];
    for i = 1:numel(all_files_data.chan)
        p = polyfit(1:numel(all_files_data.rt_all{i}),all_files_data.rt_all{i},1);
        slopes(i) = p(1);
    end
    
    f = figure();
    f.Name = 'Han_singleElectrodes_learning';
    bin_edges = [-12:2:12]*0.001;
    histogram(slopes,bin_edges);
    formatForLee(gcf)
    xlabel('Slope (RT/trial)');
    ylabel('Number of electrodes');
    set(gca,'fontsize',14);
            
            
            
            
            
            
            
            
            
            
            