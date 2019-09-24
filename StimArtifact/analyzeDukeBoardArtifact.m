%% load in a ns5


    folderpath = 'C:\Users\jts3256\Desktop\Duncan_Han_dukeProjBox\';
        
    cd(folderpath);
    file_list = dir('*.ns5');
    
    analog_pin_idx = 97;
    sync_idx = 98;
    artifact_data = {};
    sync_line_data = {};
    pwd = cd;
    window = [-4,20]; % ms
    pulse_width_1 = zeros(numel(file_list),1); % us
    pulse_width_2 = zeros(size(pulse_width_1)); % us
    interphase = 53; % us
    
    window_idx = window*30; % convert to data points

    
    for file_num = 1:numel(file_list)

        NS5 = openNSx([folderpath,file_list(file_num).name],'uV');
    
        artifact_data{file_num} = NS5.Data(analog_pin_idx,:);
        sync_line_data{file_num} = NS5.Data(sync_idx,:);
        
        % get pulse widths
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
    threshold = 0.3;
    peak_data = {};
    
    for file_num = 1:numel(file_list)
        disp(file_list(file_num).name);
        
        stim_on=find(diff(sync_line_data{file_num}-mean(sync_line_data{file_num})>3)>.5);
%         stim_on=stim_on(waveforms.waveSent == 1 & ...
%             cellfun(@isequal,waveforms.chanSent,mat2cell(analog_pin_idx+zeros(size(waveforms.chanSent)),ones(size(waveforms.chanSent)))));
%         
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
        peak_data{file_num}.cathodic.max_peak = max_cathodic_peak;
        peak_data{file_num}.cathodic.min_peak = min_cathodic_peak;

        
        [pks,locs,~,~] = findpeaks(abs(filtered_anodic_data),'MinPeakHeight',max_anodic_peak*threshold);
        locs(pks < max_anodic_peak*threshold) = [];
        locs(x_data(locs) < 0.4 | x_data(locs) > 4) = [];
        
        peak_data{file_num}.anodic.t_post_stim = x_data(locs);
        peak_data{file_num}.anodic.peak = filtered_anodic_data(locs);
        peak_data{file_num}.anodic.max_peak = max_anodic_peak;
        peak_data{file_num}.anodic.min_peak = min_anodic_peak;
        
        peak_data{file_num}.amp1 = amp_1(file_num);
        peak_data{file_num}.amp2 = amp_2(file_num);
        peak_data{file_num}.pw1 = pulse_width_1(file_num);
        peak_data{file_num}.pw2 = pulse_width_2(file_num);
        peak_data{file_num}.chan = stim_chan(file_num);
        
                %
        f=figure();
        f.Name = [file_list(file_num).name(1:end-10),'_raw'];
%         subplot(2,1,1)
        plot(x_data,(plot_data(:,1:2:6)'),'color',getColorFromList(1,0),'linewidth',1);
        hold on 
        plot(x_data,(plot_data(:,2:2:6)'),'color',getColorFromList(1,1),'linewidth',1);
        xlim([-2,6])
        ylim([-5100,5100])
        ylabel('Voltage (\muV)');
        formatForLee(gcf)
        xlabel('Time after stimulation offset (ms)');
        set(gca,'fontsize',14)
        
        saveFiguresLIB(f,folderpath,f.Name);
        
        f = figure();
%         subplot(2,1,2)
        f.Name = [file_list(file_num).name(1:end-10),'_filtered'];
        plot(x_data,(acausalFilter(plot_data(:,1:2:6))'),'color',getColorFromList(1,0),'linewidth',1);
        hold on
        plot(x_data,(acausalFilter(plot_data(:,2:2:6))'),'color',getColorFromList(1,1),'linewidth',1);
        xlim([-2,6])
        ylim([-1500,1500])
        ylabel('Voltage (\muV)');
        formatForLee(gcf)
        xlabel('Time after stimulation offset (ms)');
        set(gca,'fontsize',14)
        
        saveFiguresLIB(f,folderpath,f.Name);
        close all;
    end
    
    
%% plot amp of aux channel stim against time post stim

    colors = repmat(linspace(0,0.5,numel(peak_data))',1,3);
    
    
    for fieldname = {'cathodic','anodic'}
        f=figure(); hold on;
        f.Name = ['AmplifierRecovery_',fieldname{1}];
        for f = 1:numel(peak_data)
            ratio = peak_data{f}.(fieldname{1}).peak;
            ratio(ratio > 0) = ratio(ratio > 0)./peak_data{f}.(fieldname{1}).max_peak;
            ratio(ratio < 0) = ratio(ratio < 0)./peak_data{f}.(fieldname{1}).min_peak;
            
            
            
            plot(peak_data{f}.(fieldname{1}).t_post_stim,ratio,'color',colors(f,:))
            xlabel('Time after stimulation offset (ms)');
            ylabel('Relative gain of amplifier');
            formatForLee(gcf)
            set(gca,'fontsize',14)
        end
        
    end
    
%% make model to predict recovery as a function of time post stimulation, amplitude and polarity
    X = []; % n x m matrix
    Y = []; % n x 1 matrix

    for fieldname = {'cathodic'}
        if(strcmpi(fieldname,'cathodic'))
            is_cathodic = 1;
        else
            is_cathodic = 0;
        end
        for f = 1:numel(peak_data)
            ratio = peak_data{f}.(fieldname{1}).peak;
            ratio(ratio > 0) = ratio(ratio > 0)./peak_data{f}.(fieldname{1}).max_peak;
            ratio(ratio < 0) = ratio(ratio < 0)./peak_data{f}.(fieldname{1}).min_peak;
            t_post_stim = peak_data{f}.(fieldname{1}).t_post_stim;
            
            X = [X; t_post_stim, peak_data{f}.amp1 + zeros(size(t_post_stim)), peak_data{f}.pw1 + zeros(size(t_post_stim))];%, is_cathodic + zeros(size(t_post_stim))];
            Y = [Y;ratio'];
            
        end
        
    end
    

    % glm
    link_func = 'loglog';
    [b,dev,stats] = glmfit(X,Y,'poisson','link',link_func);
    y_fit = glmval(b,X,link_func);
    max_time_post_stim = max(t_post_stim);
    
    figure();
    plot(Y);
    hold on;
    plot(real(y_fit));
    
%% get average response and std response over all channels tested
    channels = unique(stim_chan);
    amps = unique(amp_1);
    
    gain_ratio_cathodic = cell(numel(amps),1);
    t_post_stim_cathodic = cell(numel(amps),1);
    chan_cathodic = cell(numel(amps),1);
    
    gain_ratio_anodic = cell(numel(amps),1);
    t_post_stim_anodic = cell(numel(amps),1);
    chan_anodic = cell(numel(amps),1);
    
    max_time_difference = 0.1;
    
    % for each peak data
    for p = 1:numel(peak_data)
        % find amp idx and chan idx
        amp_idx = find(amps == peak_data{p}.amp1);
        chan_idx = find(channels == peak_data{p}.chan);
        if(peak_data{p}.pw1 == 200 && peak_data{p}.pw2 == 200) % make sure this is a balanced phase case
            % store time and gain_ratio for cathodic data
            peak_vals = peak_data{p}.cathodic.peak';
            peak_vals(peak_vals > 0) = peak_vals(peak_vals > 0)/peak_data{p}.cathodic.max_peak;
            peak_vals(peak_vals < 0) = peak_vals(peak_vals < 0)/peak_data{p}.cathodic.min_peak;
            
            gain_ratio_cathodic{amp_idx} = [gain_ratio_cathodic{amp_idx};peak_vals];
            t_post_stim_cathodic{amp_idx} = [t_post_stim_cathodic{amp_idx};peak_data{p}.cathodic.t_post_stim];
            chan_cathodic{amp_idx} = [chan_cathodic{amp_idx};channels(chan_idx)+zeros(numel(peak_vals),1)];
            
            % store time and gain_ratio for anodic data
            peak_vals = peak_data{p}.anodic.peak';
            peak_vals(peak_vals > 0) = peak_vals(peak_vals > 0)/peak_data{p}.anodic.max_peak;
            peak_vals(peak_vals < 0) = peak_vals(peak_vals < 0)/peak_data{p}.anodic.min_peak;
            
            gain_ratio_anodic{amp_idx} = [gain_ratio_anodic{amp_idx};peak_vals];
            t_post_stim_anodic{amp_idx} = [t_post_stim_anodic{amp_idx};peak_data{p}.anodic.t_post_stim];
            chan_anodic{amp_idx} = [chan_anodic{amp_idx};channels(chan_idx)+zeros(numel(peak_vals),1)];
            
        end


    end

    bin_centers = 0.75:0.3:4;
    bin_width = mean(diff(bin_centers));
    min_points = 0;
    mean_gain_cathodic = nan(size(gain_ratio_cathodic,1),numel(bin_centers));
    mean_gain_anodic = nan(size(gain_ratio_anodic,1),numel(bin_centers));
    
    for a = 1:numel(t_post_stim_cathodic)
        for t = 1:numel(bin_centers)
            cathodic_keep_mask = t_post_stim_cathodic{a} > bin_centers(t)-(bin_width/2) & ...
                t_post_stim_cathodic{a} < bin_centers(t)+(bin_width/2);
            if(sum(cathodic_keep_mask) > min_points)
                mean_gain_cathodic(a,t) = mean(gain_ratio_cathodic{a}(cathodic_keep_mask),'omitnan');
            end
            anodic_keep_mask = t_post_stim_anodic{a} > bin_centers(t)-(bin_width/2) & ...
                t_post_stim_anodic{a} < bin_centers(t)+(bin_width/2);
            if(sum(anodic_keep_mask) > min_points)
                mean_gain_anodic(a,t) = mean(gain_ratio_anodic{a}(anodic_keep_mask),'omitnan');
            end
            
        end
    end
    
%     figure();
%     plot(t_post_stim_cathodic{13},gain_ratio_cathodic{13},'.','markersize',20)
%     hold on
%     plot(bin_centers,mean_gain_cathodic(13,:))
%     
    
    figure();
    ax = subplot(2,1,1);
    imagesc(bin_centers,amps,mean_gain_cathodic, [0,1]);
    ax.YDir = 'normal';
    colormap(inferno);
    formatForLee(gcf);
    ylabel('Amplitude (\muA)');
    c= colorbar(); c.Label.String = 'Gain';
    ax = subplot(2,1,2);
    imagesc(bin_centers,amps,mean_gain_anodic, [0,1]);
    ax.YDir = 'normal';
    colormap(inferno);
    formatForLee(gcf);
    xlabel('Time after stim offset (ms)')
    ylabel('Amplitude (\muA)');
    c= colorbar(); c.Label.String = 'Gain';

    
    
    
%% plot time recovery is at a certain value per amplitude (roughly)
    % interpolate mean_gain to upsample data
    bin_centers_new = bin_centers(1):0.01:bin_centers(end);
    mean_gain_cathodic_interp = interp1(bin_centers,mean_gain_cathodic',bin_centers_new)';
    mean_gain_anodic_interp = interp1(bin_centers,mean_gain_anodic',bin_centers_new)';
    
    
    gain_value = 0.9;
    amps = unique(amp_1);
    time_recovery = nan(numel(amps),2);
    
    for a = 1:numel(amps)
        t = bin_centers_new(find(mean_gain_cathodic_interp(a,:) > gain_value,1,'first'));
        if(~isempty(t))
            time_recovery(a,1) = t;
        end
        
        t = bin_centers_new(find(mean_gain_anodic_interp(a,:) > gain_value,1,'first'));
        if(~isempty(t))
            time_recovery(a,2) = t;
        end
        
    end

%     figure();
%%
    plot(amps,time_recovery(:,1),'b--','marker','.','markersize',20,'linewidth',1.5);
    hold on
    plot(amps,time_recovery(:,2),'r--','marker','.','markersize',20,'linewidth',1.5);
    
%     xlabel('Amplitude (\muA)');
%     ylabel('Amplifier recovery time after stim offset');
%     formatForLee(gcf);
%     set(gca,'fontsize',14);
%     l=legend('Cathodic','Anodic');
%     set(l,'box','off','fontsize',12,'location','best')
    
    
    