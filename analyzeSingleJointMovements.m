    input_data.folderpath = 'D:\Lab\Data\Topography\Han_20201029\';

%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Duncan_17L1\mapfiles\left S1 20190205\SN 6251-002087.cmp';
    mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Han_13B1\map files\Left S1\SN 6251-001459.cmp';
%     mapFileName = 'R:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\6250-002085\SN 6250-002085.cmp';
    
    input_data.date = '20201029';

    input_data.array = 'arrayLeftS1';
    input_data.monkey = 'monkeyHan';
    input_data.ranBy = 'ranByJoe';
    input_data.lab = 6;
    input_data.mapFile = strcat('mapFile',mapFileName);
    input_data.task = 'taskUnknown';

    pwd=cd;
    input_data.fileNum = 1;
  
    input_data.center_y = -33;
    input_data.center_x = 3;
    input_data.num_bins = 8;
    
%% make cds
    cd(input_data.folderpath)

    file_name = dir('*nev*');
    cds = commonDataStructure();
    cds.file2cds(strcat(input_data.folderpath,file_name(1).name),input_data.array,input_data.monkey,input_data.ranBy,...
        input_data.lab,input_data.mapFile,input_data.task,'recoverPreSync','ignoreJumps','ignoreFilecat');
    
    
%% import excel sheet with labeled frames, get times for each movement
    excel_sheet = dir('*frame_label*');
    frame_labels = readtable(excel_sheet(1).name);

    analog_idx = [];
    video_sync_idx = [];
    video_sync_name = 'video_sync';
    
    for i_ana = 1:numel(cds.analog)
        video_sync_idx = find(strcmpi(cds.analog{i_ana}.Properties.VariableNames,video_sync_name));
        if(~isempty(video_sync_idx))
            analog_idx = i_ana;
        end
    end    
    video_sync_data = cds.analog{analog_idx}.(video_sync_name);

    video_on = cds.analog{analog_idx}.t(find(diff(video_sync_data-mean(video_sync_data)>3)>0.5));
    video_off = cds.analog{analog_idx}.t(find(diff(video_sync_data-mean(video_sync_data)<3)>0.5));
    
    movement_id = frame_labels.Movement;
    movement_start = video_on(frame_labels.frameStart);
    movement_end = video_on(frame_labels.frameStop);
    
%% plot raster stacked by movement id for each unit
    window_ts = [-0.5,1];    
    unique_moves = uniquecell(movement_id);
 
    for i_unit = 1:numel(cds.units)
        spike_ts = cds.units(i_unit).spikes.ts;
        x_data = [];
        y_data = [];
        move_end_list = [];
        counter = 1;
        optsPlot = []; optsSave=[];
        optsPlot.DIVIDING_LINES = []; optsPlot.DIVIDING_LINES_COLORS = {};
        for i_move = 1:numel(unique_moves)
            move_idx = find(strcmpi(movement_id,unique_moves{i_move}));
            for i_trial = 1:numel(move_idx)
                % find spikes in specified channel in window
                window = movement_start(move_idx(i_trial)) + window_ts;
                spike_mask = spike_ts > window(1) & spike_ts < window(2);

                x_data = [x_data;spike_ts(spike_mask) - movement_start(move_idx(i_trial))];
                y_data = [y_data;counter*ones(sum(spike_mask),1)];
                move_end_list(end+1,1) = movement_end(move_idx(i_trial)) - movement_start(move_idx(i_trial));

                counter = counter + 1;
            end
            optsPlot.DIVIDING_LINES(end+1) = counter - 0.5;
            optsPlot.DIVIDING_LINES_COLORS{end+1} = getColorFromList(1,1);
        end
        plotRasterLIB(x_data,y_data,optsPlot,optsSave);
    end

    
%     plot(move_end_list,1:1:numel(move_end_list),'r.')

%% get firing rate during each movement for each neuron
% to look at cortical response to each movement
% baseline FR = spikes at movement start. Then find peak (Xms window)
% between movement start and end -- do (peak - baseline)/(std_baseline)
    base_window_width = [-0.2,0.2];
    peak_window_size = 0.2;
    peak_window_dt = 0.01;
    start_offset = 0.1;
    end_offset = 0.1;
    unique_moves = uniquecell(movement_id);
    
    mean_FR_move = zeros(numel(cds.units),numel(unique_moves));
    mean_FR_base = zeros(numel(cds.units),numel(unique_moves));
    std_FR_base = zeros(numel(cds.units),numel(unique_moves));
    
    row = []; col = [];
    for i_unit = 1:numel(cds.units)
        % get spikes for this unit
        spike_ts = cds.units(i_unit).spikes.ts;
        row(end+1) = cds.units(i_unit).rowNum;
        col(end+1) = cds.units(i_unit).colNum;
        for i_move = 1:numel(unique_moves)
            move_FR_trial = [];
            base_FR_trial = [];
            move_idx = find(strcmpi(movement_id,unique_moves{i_move}));
            for i_trial = 1:numel(move_idx)
                % get baseline FR
                base_window = movement_start(move_idx(i_trial)) + base_window_width;
                base_FR_trial(i_trial,1) = sum(spike_ts > base_window(1) & spike_ts < base_window(2))/diff(base_window);
                
                % get move FR
%                 spike_move_ts = spike_ts(spike_ts > movement_start(move_idx(i_trial)) & spike_ts < movement_end(move_idx(i_trial)));
%                 peak_window_start = movement_start(move_idx(i_trial));
%                 peak_windows = ((peak_window_start+start_offset):peak_window_dt:(movement_end(move_idx(i_trial))-end_offset))';
%                 peak_windows = [peak_windows, peak_windows+peak_window_size];
%                 
%                 spike_mask = spike_move_ts' > peak_windows(:,1) & spike_move_ts' < peak_windows(:,2);
%                 peak_count = max(sum(spike_mask,2));
%                 move_FR_trial(i_trial,1) = peak_count/peak_window_size;
                move_length = movement_end(move_idx(i_trial)) - movement_start(move_idx(i_trial));
                spike_mask = spike_ts > movement_start(move_idx(i_trial)) + move_length/2 - peak_window_size/2 & ...
                    spike_ts < movement_start(move_idx(i_trial)) + move_length/2 + peak_window_size/2;
                move_FR_trial(i_trial,1) = sum(spike_mask)/peak_window_size;
            end
            
            mean_FR_move(i_unit,i_move) = mean(move_FR_trial);
            mean_FR_base(i_unit,i_move) = mean(base_FR_trial);
            std_FR_base(i_unit,i_move) = std(base_FR_trial);
        end
    end
    
    z_score = (mean_FR_move - mean_FR_base)./std_FR_base;
    
    
%% plot heatmaps across the array for each movement

    num_colors = 100;
    min_z_score = -3;
    max_z_score = 3;
    base_color_pos = getColorFromList(1,0);
    base_color_neg = getColorFromList(1,1);
    
    
    colors_pos = [(0:num_colors-1)'/num_colors,(0:num_colors-1)'/num_colors,(0:num_colors-1)'/num_colors].*base_color_pos; 
    colors_neg = [(0:num_colors-1)'/num_colors,(0:num_colors-1)'/num_colors,(0:num_colors-1)'/num_colors].*base_color_neg; 
    colors = [flip(colors_neg,1);colors_pos];

    for i_moves = 1:numel(unique_moves)
        figure();
        for i_unit = 1:numel(row)            
            color_idx = min(num_colors*2,max(1,floor((z_score(i_unit,i_moves)-min_z_score)/(max_z_score-min_z_score)*(2*num_colors))));
            
            rectangle('Position',[row(i_unit),col(i_unit),1,1],'FaceColor',colors(color_idx,:), 'EdgeColor','none');
            hold on
        end
        plot([1,11,11,1,1],[1,1,11,11,1],'k','linewidth',1.5)

        f=gcf;
        set(gca,'visible','off')
        xlim([1,11])
        ylim([1,11])
        axis square

        b = colorbar;
        colormap(colors);
%         set(gca,'Visible','off');
%         b.FontSize = 14;
%         b.Ticks = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
%         b.TickDirection = 'out';
%         b.Ticks = [-1,-0.5,0,0.5,1];
%         b.Ticks = [0,0.25,0.5,0.75,1];
%         b.TickLabels = {'0','','0.5','','1'};
%         b.Label.String = 'Normalized stimulation amplitude';
%         b.Label.FontSize = 16;
%         
        
    end















    