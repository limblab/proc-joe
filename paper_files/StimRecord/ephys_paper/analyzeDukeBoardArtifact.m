
%% load in duke artifact data
load('D:\Lab\Data\stim_ephys_paper\artifact_analysis\Duncan_Han_dukeProjBox\DukeGen2_recovery_data');

%% pick a file (idx) and plot anodic and cathodic data
    make_plot = 0;
    save_plot = 0;
    threshold = 0.2;
    peak_data = {};
    num_stims = [];
    for file_num = 1:numel(file_list)
        disp(file_list(file_num).name);
        
        stim_on=find(diff(sync_line_data{file_num}-mean(sync_line_data{file_num})>3)>.5);
        
        anodic_idx = 2:2:numel(stim_on);
        cathodic_idx = 1:2:numel(stim_on);

        num_stims(end+1,:) = [numel(anodic_idx),numel(cathodic_idx)];
        
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
        
        if(make_plot && file_num==239)
            f=figure();
            f.Name = [file_list(file_num).name(1:end-10),'_raw'];
            subplot(3,1,1)
            plot(x_data,(plot_data(:,2:2:6)')/1000,'color',getColorFromList(1,1),'linewidth',1);
            hold on 
            plot(x_data,(plot_data(:,1:2:6)')/1000,'color',getColorFromList(1,0),'linewidth',0.5);
            xlim([-2,5])
            ylim([-5.1,5.1])
            ylabel('Voltage (mV)');
            formatForLee(gcf)
            set(gca,'fontsize',14)

            subplot(3,1,2)
            f.Name = [file_list(file_num).name(1:end-10),'_filtered'];
            plot(x_data,(acausalFilter(plot_data(:,2:2:6))')/1000,'color',getColorFromList(1,1),'linewidth',1);
            hold on
            plot(x_data,(acausalFilter(plot_data(:,1:2:6))')/1000,'color',getColorFromList(1,0),'linewidth',0.5);
            xlim([-2,5])
            ylim([-2,2])
            ylabel('Voltage (mV)');
            formatForLee(gcf)
            set(gca,'fontsize',14)
            
            
            
            % this assumes that the rest of the script has been ran, which is dumb but it's a hack...
            gain_data_cath = ones(size(x_data));
            gain_data_ano = ones(size(x_data));
            amp_idx = find(amps==50);
            t_data = x_data;
            
            gain_data_cath_interp = interp1(bin_centers-0.05, metric_gain_cathodic(amp_idx,:),t_data(t_data > 0 & t_data < 5),'linear','extrap');
            gain_data_ano_interp = interp1(bin_centers-0.05, metric_gain_anodic(amp_idx,:),t_data(t_data > 0 & t_data < 5),'linear','extrap');
            
            gain_data_cath_interp(gain_data_cath_interp < 0) = 0.00001;
            gain_data_ano_interp(gain_data_ano_interp < 0) = 0.00001;
            
            gain_data_cath_interp(gain_data_cath_interp > 1) = 1;
            gain_data_ano_interp(gain_data_ano_interp > 1) = 1;
            
            gain_data_cath(t_data>0 & t_data < 5) = gain_data_cath_interp;
            gain_data_ano(t_data>0 & t_data < 5) = gain_data_ano_interp;
           
            subplot(3,1,3)
            f.Name = [file_list(file_num).name(1:end-10),'_filtered'];
            % anodic
            plot(x_data,1./(gain_data_cath)'.*(acausalFilter(plot_data(:,2:2:6))')/1000,'color',getColorFromList(1,1),'linewidth',1);
            hold on
            % cathodic
            plot(x_data,1./(gain_data_ano)'.*(acausalFilter(plot_data(:,1:2:6))')/1000,'color',getColorFromList(1,0),'linewidth',0.5);
            xlim([-2,5])
            ylim([-2,2])
            ylabel('Voltage (mV)');
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
    bin_centers = 0.75:0.5:4;
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
    
    amp_keep = [10,25,50,100]';
    amp_mask = any(amps==amp_keep);
    
    offset_cathodic = [-0.8,0.0,1.2,-0.8]./40;
    offset_anodic = [-1.2,-0.2,-0.8,-1.2]./40;
    color_list = inferno(numel(amp_keep)+1);
    
    f=figure(); hold on;
    f.Name = 'AmplifierRecovery_gainVsTime';
    alpha_min = 0.4;
    
    for i_amp = 1:numel(amps)
        for polarity_idx = 1:2
            if(amp_mask(i_amp)==1)
                if(polarity_idx == 1)
%                     low_bar = metric_gain_cathodic(i_amp,:)-squeeze(prc_gain_cathodic(i_amp,:,1));
%                     high_bar = squeeze(prc_gain_cathodic(i_amp,:,2))-metric_gain_cathodic(i_amp,:);
                    
                    low_bar = -std_gain_cathodic(i_amp,:);
                    high_bar = std_gain_cathodic(i_amp,:);
%                     plot(bin_centers,metric_gain_cathodic(i_amp,:),...
                    errorbar(bin_centers+offset_cathodic(sum(amp_mask(1:i_amp))),metric_gain_cathodic(i_amp,:),low_bar,high_bar,...
                        'marker','.','markersize',22,'linestyle','-',...
                        'color',color_list(sum(amp_mask(1:i_amp)),:),'linewidth',2);
                    
%                     
        

                else
%                     low_bar = metric_gain_anodic(i_amp,:)-squeeze(prc_gain_anodic(i_amp,:,1));
%                     high_bar = squeeze(prc_gain_anodic(i_amp,:,2))-metric_gain_anodic(i_amp,:);
                    low_bar = -std_gain_anodic(i_amp,:);
                    high_bar = std_gain_anodic(i_amp,:);
                    d=errorbar(bin_centers+offset_anodic(sum(amp_mask(1:i_amp))),metric_gain_anodic(i_amp,:),low_bar,high_bar,...
                        'marker','s','markersize',8,'linestyle','--',...
                        'color',color_list(sum(amp_mask(1:i_amp)),:),'linewidth',2);
                    d.Bar.LineStyle = 'dotted';
                %                     plot(bin_centers,metric_gain_anodic(i_amp,:),...

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

    xlim([0,5]);
    ylim([0.2,1.15]);
%% stats on gain data -- with gain data with a curve, let pulse polarity and amplitude be factors

% amps has the amplitude, gain_ratio has the cathodic and anodic data for
% each time point, t_post_stim has the time points for that data

% format data into a table. 

t_list = []; gain_list = []; is_cathodic = []; amp_list = []; stim_chan_list = [];
for i_amp = 1:numel(gain_ratio_anodic)
    for i_case = 1:2
        switch i_case
            case 1
                temp_gain = gain_ratio_anodic{i_amp};
                temp_t = t_post_stim_anodic{i_amp};
                temp_chan = chan_anodic{i_amp};
                temp_is_cath = 0;
            case 2
                temp_gain = gain_ratio_cathodic{i_amp};
                temp_t = t_post_stim_cathodic{i_amp};
                temp_chan = chan_cathodic{i_amp};
                temp_is_cath = 1;
        end
        
        keep_mask = temp_t > 0.9 & temp_t < 1.1;
%         keep_mask = temp_t > 1.3 & temp_t < 1.4;
%         keep_mask = temp_t > 2.0 & temp_t < 2.1;
%         keep_mask = ones(size(temp_t));
        t_list = [t_list; temp_t(keep_mask==1)];
        gain_list = [gain_list; temp_gain(keep_mask==1)];
        amp_list = [amp_list; amps(i_amp)*ones(sum(keep_mask==1),1)];
        is_cathodic = [is_cathodic; temp_is_cath*ones(sum(keep_mask==1),1)];
        stim_chan_list = [stim_chan_list; temp_chan(keep_mask==1)];
    end
    
end

gain_tbl = table(t_list,gain_list,amp_list,is_cathodic,stim_chan_list,'VariableNames',{'t','gain','amp','is_cath','chan'});

gain_tbl.is_cath = categorical(gain_tbl.is_cath);
gain_tbl.chan = categorical(gain_tbl.chan);


mdl_spec = 'gain~amp+is_cath + chan';

gain_mdl = fitlm(gain_tbl,mdl_spec)
% gain_mdl = fitglm(gain_tbl,mdl_spec,'Distribution','binomial')
gain_pred = predict(gain_mdl,gain_tbl);
plot(gain_pred,gain_tbl.gain,'.')


