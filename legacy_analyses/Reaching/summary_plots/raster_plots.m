%% load in cds

load('D:\Lab\Data\FreeReaching\Han_20210623\neural-data\Han_20210623_cds.mat');
% load('D:\Lab\Data\FreeReaching\Rocket_20210723\neural-data\Rocket_20210723_cds.mat');

%% pick random times within dataWindow which are not overlapping

num_trials = 6;
time_plot = 4;

optsPlot = []; optsSave = [];
optsPlot.MARKER_STYLE = 'line'; % line or .
optsPlot.X_LIMITS = [0,time_plot];
optsPlot.LINE_WIDTH = 0.5;
optsPlot.MAKE_FIGURE = 0;

ana_idx = 3;


for i_cds = 1:2
    hand_vel_data = [];
    hand_pos_data = [cds_list{i_cds}.analog{ana_idx}.hand2_x,cds_list{i_cds}.analog{ana_idx}.hand2_y,cds_list{i_cds}.analog{ana_idx}.hand2_z];
    t_data = cds_list{i_cds}.analog{ana_idx}.t;
    for i=1:3
        hand_vel_data(:,i) = gradient(hand_pos_data(:,i), mean(diff(t_data)));
        hand_vel_data(:,i) = fillmissing(hand_vel_data(:,i),'linear');
    end
    hand_spd_data = sqrt(sum(hand_vel_data.^2,2));
    
    times_started = [];
    for i_trial = 1:num_trials
        % get time to start plotting
        is_overlapping = 1;
        while(is_overlapping)
            time_start = rand()*(diff(cds_list{i_cds}.meta.dataWindow)-time_plot) + cds_list{i_cds}.meta.dataWindow(1);
            if(i_trial > 1)
                is_overlapping = any(times_started - time_start <= time_plot & times_started - time_start > 0);
            else 
                is_overlapping = 0;
            end
        end
        times_started(end+1) = time_start;

        % make raster plot
        unit_plot = 1;
        x_data = [];
        y_data = [];

        f=figure('Position',[2283 550 220 250]);
        f.Name = ['Han_20210623_rasterRaster_trial', num2str(i_trial),'_task',cds_list{i_cds}.meta.task];
        for i_unit = 1:numel(cds_list{i_cds}.units)
            if(cds_list{i_cds}.units(i_unit).ID >= 1)
                keep_mask = cds_list{i_cds}.units(i_unit).spikes.ts >= time_start & cds_list{i_cds}.units(i_unit).spikes.ts <= time_start + time_plot;
                x_data = [x_data;cds_list{i_cds}.units(i_unit).spikes.ts(keep_mask) - time_start];
                y_data = [y_data; unit_plot*ones(sum(keep_mask),1)];


                unit_plot = unit_plot + 1;
            end
        end
        output_data = plotRasterLIB(x_data,y_data,optsPlot,optsSave);

        % plot hand speed during trial
        f=figure('Position',[2283 115 220 250]);
        f.Name = ['Han_20210623_speedRaster_trial', num2str(i_trial),'_task',cds_list{i_cds}.meta.task];
        hand_vel_idx = [find(t_data >= time_start,1,'first'), find(t_data >= time_start+time_plot,1,'first')];
        plot(linspace(0,time_plot,diff(hand_vel_idx)+1),hand_spd_data(hand_vel_idx(1):hand_vel_idx(2)));
        formatForLee(gcf);
    end
end