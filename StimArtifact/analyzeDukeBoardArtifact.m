%% load in a ns5
    
    folderpath = 'E:\Data\Joseph\Duncan_stim_data\Duncan_20191124_longTrains_dukeRecovery\';
        
    cd(folderpath);
    file_list = dir('*.ns5');
    
    analog_pin_idx = 97;
    sync_idx = 98;
    artifact_data = {};
    sync_line_data = {};
    pwd = cd;
    window = [-4,5]; % ms
    pulse_width_1 = zeros(numel(file_list),1); % us
    pulse_width_2 = zeros(size(pulse_width_1)); % us
    interphase = 53; % us
    
    window_idx = window*30; % convert to data points
    
    for file_num = 1:numel(file_list)
        NS5 = openNSx([folderpath,file_list(file_num).name],'uV');
%         if(size(NS5.Data,2) > 500000)
%             data = NS5.Data(1,1:500000);
%             save([file_list(file_num).name(1:end-4),'_data'],'data');
%         end
%         f=figure();
%         f.Name = file_list(file_num).name;
%         periodogram(NS5.Data(1,:))
%         hold on
        artifact_data{file_num} = NS5.Data(analog_pin_idx,:);
        sync_line_data{file_num} = NS5.Data(sync_idx,:);
%         
%         % get pulse widths
        pw1_idx = strfind(file_list(file_num).name,'PW1');
        pw2_idx = strfind(file_list(file_num).name,'PW2');
        amp1_idx = strfind(file_list(file_num).name,'A1');
        amp2_idx = strfind(file_list(file_num).name,'A2');
        stim_idx = strfind(file_list(file_num).name,'stim');
        chan_idx = strfind(file_list(file_num).name,'chan');
        underscore_idx = strfind(file_list(file_num).name,'_');
        
        pulse_width_1(file_num) = str2num(file_list(file_num).name(pw1_idx+4:underscore_idx(find(underscore_idx > pw1_idx,1,'first'))-1));
        pulse_width_2(file_num) = str2num(file_list(file_num).name(pw2_idx+4:underscore_idx(find(underscore_idx > pw2_idx,1,'first'))-1));
        amp_1(file_num) = str2num(file_list(file_num).name(amp1_idx+3:underscore_idx(find(underscore_idx > amp1_idx,1,'first'))-1));
        amp_2(file_num) = str2num(file_list(file_num).name(amp2_idx+3:underscore_idx(find(underscore_idx > amp2_idx,1,'first'))-1));
        stim_chan(file_num) = str2num(file_list(file_num).name(chan_idx+4:stim_idx-1));
    end
    cd(pwd);
    
%% pick a file (idx) and plot anodic and cathodic data
    make_plot = 0;
    save_plot = 0;
    threshold = 0.2;
    peak_data = {};
    
    for file_num = 1:numel(file_list)
        disp(file_list(file_num).name);
        
        stim_on=find(diff(sync_line_data{file_num}-mean(sync_line_data{file_num})>3)>.5);

        anodic_idx = 2:2:numel(stim_on);
        cathodic_idx = 1:2:numel(stim_on);

        plot_data = [];
        x_data = [window_idx(1):window_idx(2)]'/30 - (pulse_width_1(file_num) + pulse_width_2(file_num) + interphase)/1000;

        for st = 1:numel(stim_on)
            plot_data(:,st) = artifact_data{file_num}((stim_on(st)+window_idx(1)):(stim_on(st)+window_idx(2)));
        end
        
        filtered_cathodic_data = mean(acausalFilter(plot_data(:,cathodic_idx))');
        filtered_anodic_data = mean(acausalFilter(plot_data(:,anodic_idx))');
        
        max_cathodic_peak = max(filtered_cathodic_data(find(x_data > 4,1,'first'):find(x_data > 7,1,'first')));
        min_cathodic_peak = min(filtered_cathodic_data(find(x_data > 4,1,'first'):find(x_data > 7,1,'first')));
        
        max_anodic_peak = max(filtered_anodic_data(find(x_data > 4,1,'first'):find(x_data > 7,1,'first')));
        min_anodic_peak = min(filtered_cathodic_data(find(x_data > 4,1,'first'):find(x_data > 7,1,'first')));
        
        [pks,locs,~,~] = findpeaks(abs(filtered_cathodic_data),'MinPeakHeight',max_cathodic_peak*threshold);
        locs(pks < max_cathodic_peak*threshold) = [];
        locs(x_data(locs) < 0.4 | x_data(locs) > 4) = [];
        
        peak_data{file_num}.cathodic.t_post_stim = x_data(locs);
        peak_data{file_num}.cathodic.peak = filtered_cathodic_data(locs);
        peak_data{file_num}.cathodic.max_peak = mean(peak_data{file_num}.cathodic.peak(peak_data{file_num}.cathodic.t_post_stim' > 3 & peak_data{file_num}.cathodic.peak > 0));
        peak_data{file_num}.cathodic.min_peak = mean(peak_data{file_num}.cathodic.peak(peak_data{file_num}.cathodic.t_post_stim' > 3 & peak_data{file_num}.cathodic.peak < 0));

        [pks,locs,~,~] = findpeaks(abs(filtered_anodic_data),'MinPeakHeight',max_anodic_peak*threshold);
        locs(pks < max_anodic_peak*threshold) = [];
        locs(x_data(locs) < 0.4 | x_data(locs) > 4) = [];
        
        peak_data{file_num}.anodic.t_post_stim = x_data(locs);
        peak_data{file_num}.anodic.peak = filtered_anodic_data(locs);
        peak_data{file_num}.anodic.max_peak = mean(peak_data{file_num}.anodic.peak(peak_data{file_num}.anodic.t_post_stim' > 3 & peak_data{file_num}.anodic.peak > 0));
        peak_data{file_num}.anodic.min_peak = mean(peak_data{file_num}.anodic.peak(peak_data{file_num}.anodic.t_post_stim' > 3 & peak_data{file_num}.anodic.peak < 0));
        
        peak_data{file_num}.amp1 = amp_1(file_num);
        peak_data{file_num}.amp2 = amp_2(file_num);
        peak_data{file_num}.pw1 = pulse_width_1(file_num);
        peak_data{file_num}.pw2 = pulse_width_2(file_num);
        peak_data{file_num}.chan = stim_chan(file_num);
        
        if(make_plot)
            f=figure();
            f.Name = [file_list(file_num).name(1:end-10),'_raw'];
            subplot(2,1,1)
            plot(x_data,(plot_data(:,1:2:6)'),'color',getColorFromList(1,0),'linewidth',1);
            hold on 
            plot(x_data,(plot_data(:,2:2:6)'),'color',getColorFromList(1,1),'linewidth',1);
            xlim([-2,6])
            ylim([-5100,5100])
            ylabel('Voltage (\muV)');
            formatForLee(gcf)
            xlabel('Time after stimulation offset (ms)');
            set(gca,'fontsize',14)

            subplot(2,1,2)
            f.Name = [file_list(file_num).name(1:end-10),'_filtered'];
            plot(x_data,(acausalFilter(plot_data(:,1:2:6))'),'color',getColorFromList(1,0),'linewidth',1);
            hold on
            plot(x_data,(acausalFilter(plot_data(:,2:2:6))'),'color',getColorFromList(1,1),'linewidth',1);
            xlim([-2,6])
            ylim([-3000,3000])
            ylabel('Voltage (\muV)');
            formatForLee(gcf)
            xlabel('Time after stimulation offset (ms)');
            set(gca,'fontsize',14)
            
            if(save_plot)
                saveFiguresLIB(f,folderpath,f.Name);
                close all;
            end
        end
    end
    
    

    
%% get average response and std response over all channels tested
    channels = unique(stim_chan);
    amps = unique(amp_1);
    gain_value = 0.75;
    use_peak_to_peak = 1;
    use_only_pos_vals = 0;
    
    gain_ratio_cathodic = cell(numel(amps),1);
    t_post_stim_cathodic = cell(numel(amps),1);
    chan_cathodic = cell(numel(amps),1);
    
    gain_ratio_anodic = cell(numel(amps),1);
    t_post_stim_anodic = cell(numel(amps),1);
    chan_anodic = cell(numel(amps),1);
    
    max_time_difference = 0.15;
    
    
    interp_times = 0.5:0.01:4;
    
    time_recovered = cell(size(gain_ratio_cathodic,1),2); % cathodic, anodic
    
    % for each peak data (electrode stimulated at a certain amplitude)
    for p = 1:numel(peak_data)
        % find amp idx and chan idx
        amp_idx = find(amps == peak_data{p}.amp1);
        chan_idx = find(channels == peak_data{p}.chan);
        if(peak_data{p}.pw1 == 200 && peak_data{p}.pw2 == 200) % make sure this is a balanced phase case
            % store time and gain_ratio for cathodic data
            if(use_peak_to_peak)
                peak_vals_cathodic = []; peak_vals_anodic = [];
                t_cathodic = []; t_anodic = [];
                
                for i_peak = 1:numel(peak_data{p}.cathodic.peak)
                    if(peak_data{p}.cathodic.peak(i_peak) > 0) % try to find closest min peak
                        t_list = peak_data{p}.cathodic.t_post_stim;
                        t_list(i_peak) = Inf;
                        [dt,closest_min_peak] = min(abs(t_list - peak_data{p}.cathodic.t_post_stim(i_peak)));
                        if(dt < max_time_difference)
                            peak_vals_cathodic(end+1,1) = peak_data{p}.cathodic.peak(i_peak) - peak_data{p}.cathodic.peak(closest_min_peak);
                            t_cathodic(end+1,1) = peak_data{p}.cathodic.t_post_stim(i_peak) + dt/2;
                        end
                    end
                end
                peak_vals_cathodic = peak_vals_cathodic/(peak_data{p}.cathodic.max_peak - peak_data{p}.cathodic.min_peak);
                
                for i_peak = 1:numel(peak_data{p}.anodic.peak)
                    if(peak_data{p}.anodic.peak(i_peak) > 0) % try to find closest min peak
                        t_list = peak_data{p}.anodic.t_post_stim;
                        t_list(i_peak) = Inf;
                        [dt,closest_min_peak] = min(abs(t_list - peak_data{p}.anodic.t_post_stim(i_peak)));
                        if(dt < max_time_difference)
                            peak_vals_anodic(end+1,1) = peak_data{p}.anodic.peak(i_peak) - peak_data{p}.anodic.peak(closest_min_peak);
                            t_anodic(end+1,1) = peak_data{p}.anodic.t_post_stim(i_peak) + dt/2;
                        end
                    end
                end
                peak_vals_anodic = peak_vals_anodic/(peak_data{p}.anodic.max_peak - peak_data{p}.anodic.min_peak);
                
            elseif(use_only_pos_vals) % look at only positive peaks
                peak_vals_cathodic = peak_data{p}.cathodic.peak(peak_data{p}.cathodic.peak>0)'/peak_data{p}.cathodic.max_peak;
                peak_vals_anodic = peak_data{p}.anodic.peak(peak_data{p}.anodic.peak>0)'/peak_data{p}.anodic.max_peak;
                t_cathodic = peak_data{p}.cathodic.t_post_stim(peak_data{p}.cathodic.peak>0);
                t_anodic = peak_data{p}.anodic.t_post_stim(peak_data{p}.anodic.peak>0);
                
            else % look at positive and negative peaks separately
                peak_vals_cathodic = peak_data{p}.cathodic.peak';
                peak_vals_cathodic(peak_vals_cathodic > 0) = peak_vals_cathodic(peak_vals_cathodic > 0)/peak_data{p}.cathodic.max_peak;
                peak_vals_cathodic(peak_vals_cathodic < 0) = peak_vals_cathodic(peak_vals_cathodic < 0)/peak_data{p}.cathodic.min_peak;
                
                peak_vals_anodic = peak_data{p}.anodic.peak';
                peak_vals_anodic(peak_vals_anodic > 0) = peak_vals_anodic(peak_vals_anodic > 0)/peak_data{p}.anodic.max_peak;
                peak_vals_anodic(peak_vals_anodic < 0) = peak_vals_anodic(peak_vals_anodic < 0)/peak_data{p}.anodic.min_peak;
                
                t_cathodic = peak_data{p}.cathodic.t_post_stim;
                t_anodic = peak_data{p}.anodic.t_post_stim;
            end
            
            % store time and gain_ratio for cathodic and anodic data
            gain_ratio_cathodic{amp_idx} = [gain_ratio_cathodic{amp_idx};peak_vals_cathodic];
            t_post_stim_cathodic{amp_idx} = [t_post_stim_cathodic{amp_idx};t_cathodic];
            chan_cathodic{amp_idx} = [chan_cathodic{amp_idx};channels(chan_idx)+zeros(numel(peak_vals_cathodic),1)];
            
            gain_ratio_anodic{amp_idx} = [gain_ratio_anodic{amp_idx};peak_vals_anodic];
            t_post_stim_anodic{amp_idx} = [t_post_stim_anodic{amp_idx};t_anodic];
            chan_anodic{amp_idx} = [chan_anodic{amp_idx};channels(chan_idx)+zeros(numel(peak_vals_anodic),1)];
            
            % get time recovered for each channel
            recovery_cathodic = interp1(t_cathodic,peak_vals_cathodic,interp_times)';
            recovery_anodic = interp1(t_anodic,peak_vals_anodic,interp_times)';
            
            time_recovered{amp_idx,1}(end+1,1) = interp_times(find(recovery_cathodic > gain_value,1,'first'));
            time_recovered{amp_idx,2}(end+1,1) = interp_times(find(recovery_anodic > gain_value,1,'first'));
        end


    end

%% plot gain of amplifier at different time across amplitudes and polarity
    bin_centers = 0.75:0.5:3.75;
    bin_width = mean(diff(bin_centers));
    min_points = 0;
    metric_gain_cathodic = nan(size(gain_ratio_cathodic,1),numel(bin_centers));
    metric_gain_anodic = nan(size(gain_ratio_anodic,1),numel(bin_centers));
    
    std_gain_cathodic = nan(size(metric_gain_cathodic));
    std_gain_anodic = nan(size(metric_gain_anodic));
    prc_gain_cathodic = nan(size(metric_gain_cathodic,1),size(metric_gain_cathodic,2),2);
    prc_gain_anodic = nan(size(metric_gain_anodic,1),size(metric_gain_anodic,2),2);
    for a = 1:numel(t_post_stim_cathodic)
        for t = 1:numel(bin_centers)
            cathodic_keep_mask = t_post_stim_cathodic{a} > bin_centers(t)-(bin_width/2) & ...
                t_post_stim_cathodic{a} < bin_centers(t)+(bin_width/2);
            if(sum(cathodic_keep_mask) > min_points)
                metric_gain_cathodic(a,t) = mean(gain_ratio_cathodic{a}(cathodic_keep_mask),'omitnan');
                std_gain_cathodic(a,t) = std(gain_ratio_cathodic{a}(cathodic_keep_mask),'omitnan');
                prc_gain_cathodic(a,t,:) = prctile(gain_ratio_cathodic{a}(cathodic_keep_mask),[5,95]);
            end
            anodic_keep_mask = t_post_stim_anodic{a} > bin_centers(t)-(bin_width/2) & ...
                t_post_stim_anodic{a} < bin_centers(t)+(bin_width/2);
            if(sum(anodic_keep_mask) > min_points)
                metric_gain_anodic(a,t) = mean(gain_ratio_anodic{a}(anodic_keep_mask),'omitnan');
                std_gain_anodic(a,t) = std(gain_ratio_anodic{a}(anodic_keep_mask),'omitnan');
                prc_gain_anodic(a,t,:) = prctile(gain_ratio_anodic{a}(anodic_keep_mask),[5,95]);
            end
            
        end
    end
    
    amp_keep = [15,30,60,100]';
    offset = (linspace(0,1,numel(amp_keep)))*0.05;
    amp_mask = any(amps==amp_keep);
    
    color_list = inferno(numel(amp_keep)+1);
    
    f=figure(); hold on;
    f.Name = 'AmplifierRecovery_gainVsTime';
    alpha_min = 0.4;
    
    for i_amp = 1:numel(amps)
        for polarity_idx = 1:2
            if(amp_mask(i_amp)==1)
                if(polarity_idx == 1)
%                     low_bar = median_gain_cathodic(i_amp,:)-squeeze(prc_gain_cathodic(i_amp,:,1));
%                     high_bar = squeeze(prc_gain_cathodic(i_amp,:,2))-median_gain_cathodic(i_amp,:);
                    
%                     errorbar(bin_centers+offset(sum(amp_mask(1:i_amp))),median_gain_cathodic(i_amp,:),low_bar,high_bar,...
                    plot(bin_centers,metric_gain_cathodic(i_amp,:),...
                        'marker','.','markersize',12,'linestyle','-',...
                        'color',color_list(sum(amp_mask(1:i_amp)),:),'linewidth',2);
                else
%                     low_bar = median_gain_cathodic(i_amp,:)-squeeze(prc_gain_cathodic(i_amp,:,1));
%                     high_bar = squeeze(prc_gain_cathodic(i_amp,:,2))-median_gain_cathodic(i_amp,:);
%                     errorbar(bin_centers-offset(sum(amp_mask(1:i_amp))),median_gain_anodic(i_amp,:),low_bar,high_bar,...
                    plot(bin_centers,metric_gain_anodic(i_amp,:),...
                        'marker','.','markersize',12,'linestyle','--',...
                        'color',color_list(sum(amp_mask(1:i_amp)),:),'linewidth',2);
                end
            end
        end
    end
    ax = gca;
    ax.XLim(1) = 0;    
    formatForLee(gcf);
    set(ax,'fontsize',14);
    xlabel('Time after stim offset (ms)');
    ylabel('Relative gain');
    
    l=legend('Cathodic-first','Anodic-first');
    set(l,'box','off');

%% plot time recovery is at a certain value per amplitude
    % interpolate mean_gain to upsample data
    mean_time_recovery = cellfun(@mean,time_recovered);
    std_time_recovery = cellfun(@std,time_recovered);
    amp_offset = [-0.5,0.5];
    
    amps = unique(amp_1);
    
    figure();
    errorbar(amps+amp_offset(1),mean_time_recovery(:,1),std_time_recovery(:,1),...
        '--','color',getColorFromList(1,0),'marker','.','markersize',20,'linewidth',1.5);
    hold on
    errorbar(amps+amp_offset(2),mean_time_recovery(:,2),std_time_recovery(:,2),...
        '--','color',getColorFromList(1,1),'marker','.','markersize',20,'linewidth',1.5);
    
    xlabel('Amplitude (\muA)');
    ylabel('Amplifier recovery time after stim offset');
    formatForLee(gcf);
    set(gca,'fontsize',14);
    l=legend('Cathodic','Anodic');
    set(l,'box','off','fontsize',12,'location','best')
    
    
% %% plot heatmap of gain recovery for anodic vs cathodic
% % make plot depicting gain recovery over time across amplitudes
%     for a = 1:numel(t_post_stim_cathodic)
%         for t = 1:numel(bin_centers)
%             cathodic_keep_mask = t_post_stim_cathodic{a} > bin_centers(t)-(bin_width/2) & ...
%                 t_post_stim_cathodic{a} < bin_centers(t)+(bin_width/2);
%             if(sum(cathodic_keep_mask) > min_points)
%                 metric_gain_cathodic(a,t) = median(gain_ratio_cathodic{a}(cathodic_keep_mask),'omitnan');
%             end
%             anodic_keep_mask = t_post_stim_anodic{a} > bin_centers(t)-(bin_width/2) & ...
%                 t_post_stim_anodic{a} < bin_centers(t)+(bin_width/2);
%             if(sum(anodic_keep_mask) > min_points)
%                 metric_gain_anodic(a,t) = median(gain_ratio_anodic{a}(anodic_keep_mask),'omitnan');
%             end
%             
%         end
%     end
%     
%     amp_step = 10;
%     amps_mask = mod(amps,amp_step) == 0;
%     
%     figure();
%     for polarity_idx = 1:2
%         ax = subplot(2,1,polarity_idx);
%         if(polarity_idx == 1)
%             imagesc(bin_centers,amps(amps_mask),metric_gain_cathodic(amps_mask,:), [0,1]);
%         else
%             imagesc(bin_centers,amps(amps_mask),metric_gain_anodic(amps_mask,:), [0,1]);
%             xlabel('Time after stim offset (ms)')
%         end
%         
%         ax.YDir = 'normal';
%         ax.XTick = bin_centers(2:2:end);
%         ax.XAxis.MinorTickValues = bin_centers;
%     
%         ax.YAxis.MinorTickValues = amps(amps_mask);
%         ax.YTick = ax.YAxis.MinorTickValues(1:3:end);
%         
%         colormap(inferno);
%         formatForLee(gcf);
%         ylabel('Amplitude (\muA)');
%         c= colorbar(); c.Label.String = 'Gain';
%         c.Ticks = [0:0.25:1];
%         c.TickLength = 0.03;
%         c.TickDirection = 'out';
%         for tick_idx = 2:2:numel(c.Ticks)
%             c.TickLabels{tick_idx} = '';
%         end
%         c.FontSize = 14;
%         set(ax,'fontsize',14)
%         
%     end